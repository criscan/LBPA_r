#' @title Length-Based Pseudo-cohort Analysis - LBPA
#' @description The model estimates the fishing mortality and depletion of the Spawning Potential Ratio (SPR) in an equilibrium condition.
#' The model is fitted to several length frequencies of catches simultaneously by minimizing a penalized log-likelihood function.
#' @param name excel (.xlsx) file with data and parameters, graph_opt (=T or F) is graphics display is required, and save_opt (=T or F) if the csv files are generated: model parameters,
#' population variables, log-likelihood and per-recruit variables . The xlsx file contains four sheets: LF data, biological parameters,
#' initial parameters and data weighting
#' @return Several variables with model fit and population variables. 
#' @examples  lbpaout=LBPA_fits("name.xlsx", graph_opt=T, save_opt=T)


LBPA_fits <- function(name, graph_opt, save_opt) {
  library(readxl)
  library(openxlsx)
  library(areaplot)
  
  # Lectura de datos desde archivo Excel
  datos    <- as.matrix(read_xlsx(name, sheet = 1, col_names = TRUE))
  parbiol  <- as.numeric(read_xlsx(name, sheet = 2, col_names = TRUE))
  pars_ini <- as.numeric(read_xlsx(name, sheet = 3, col_names = TRUE)[1,])
  cvpar    <- as.numeric(read_xlsx(name, sheet = 3, col_names = TRUE)[2,])
  lambda   <- as.numeric(read_xlsx(name, sheet = 4, col_names = TRUE))
  
  if(length(parbiol)<10){stop("inssuficient biological parameters")}
  
  
  # Parámetros biológicos y valores iniciales
  Talla <- datos[, 1]
  M     <- parbiol[3]
  L50   <- pars_ini[1]; slope <- pars_ini[2]; Fcr <- pars_ini[3]
  Lr    <- pars_ini[4]; a0    <- pars_ini[5];  cv  <- pars_ini[6]
  target <- parbiol[10]
  parini <- log(c(L50, slope, Fcr, Lr, a0 + 1e-5, cv + 1e-5))
  
  # Tamaño de muestra y lista de entrada
  nm <- 500
  data_list <- list(
    datos = datos,
    parbiol = parbiol,
    nm = nm,
    cv_par = cvpar,
    prioris = parini,
    lambda = lambda,
    flag=1)
  
  # Función de verosimilitud-----------------------------------------------------------
  ferrLBPA <- function(parini, data_list) {
    with(as.list(data_list), {
      
      par_vals <- exp(parini)
      
      names(par_vals) <- c("L50", "slope", "Fcr", "Lr", "a0", "cv")
      
      Tmax <- round(-log(0.01)/parbiol[3])
      
      key <- exp(prioris[4]) / parbiol[1]
         tr  <- round(-1 / parbiol[2] * log(1 - exp(prioris[4]) / parbiol[1]) - 0.5)
    Lmm=round(2*parbiol[6]-parbiol[7]) #length at 5% maturity
    tmm  <- round(-1 / parbiol[2] * log(1 - Lmm / parbiol[1]) - 0.5)
    target<-parbiol[10]
    
    age=c(min(tmm,tr):Tmax) #select the minimum between tr and tmm
    Nages=length(age)


    Lage <- numeric(Nages)
    Lage[1] <- par_vals["Lr"]
    for (i in 2:Nages) {
      Lage[i] <- Lage[i - 1] * exp(-parbiol[2]) + parbiol[1] * (1 - exp(-parbiol[2]))
    }
    
    sd_age <- par_vals["a0"] + par_vals["cv"] * Lage
    Sel    <- 1 / (1 + exp(-log(19) * (Lage - par_vals["L50"]) / par_vals["slope"]))
    Fmort <- par_vals["Fcr"] * Sel
    Z     <- Fmort + parbiol[3]
    
  
    N  <- rep(1, Nages)
    for (i in 2:Nages) {
      N[i]  <- N[i - 1] * exp(-Z[i - 1])
    }
    N[Nages]  <- N[Nages] / (1 - exp(-Z[Nages]))
    
    C  <- N * Fmort / Z * (1 - exp(-Z))
    
    pdf=pdf_mat(Talla,Lage, sd_age)
    
    Ctalla <- t(pdf) %*% C
    ppred  <- Ctalla / sum(Ctalla)
    
    fobs <- datos[, -1]
    lsum <- 0
    for (i in 2:ncol(datos)) {
      lsum <- lsum + nm * lambda[i - 1] * sum(fobs[, i - 1] / sum(fobs[, i - 1]) * log(ppred+1e-10))
    }
    
    priors=(prioris - log(par_vals))^2 / cv_par^2
    penal <- 0.5 * sum(priors)
    fun<--lsum + penal
    
    if (flag==1) {
      return(fun)
    }
    
    if(flag==2){
      pobs=fobs
      
      for (i in 1:ncol(pobs)){
        pobs[,i]=pobs[,i]/sum(pobs[,i])
      }
      
    }
    
    other_1=list(lprop=fun,
                 lprior=priors,
                 talla=datos[, 1],
                 ppred=ppred,
                 pobs=pobs,
                 Lage=Lage,
                 sd_age=sd_age,
                 ages=age,
                 Nages=Nages,
                 Tmax=Tmax,
                 Lr=Lr,
                 tr=tr,
                 Lmm=Lmm,
                 tmm=tmm,
                 tar=target)

    return(other_1)
    
  })}

# Funcion pdf---------------

pdf_mat<-function(Talla,Lage, sd_age){

# Matriz talla-edad inversa
Nages=length(Lage)
dtalla <- Talla[2] - Talla[1]
xs     <- 0.5 * dtalla
pdf    <- matrix(NA, Nages, length(Talla))
for (i in 1:Nages) {
  for (j in 1:length(Talla)) {
    z1 <- ((Talla[j] - xs) - Lage[i]) / sd_age[i]
    z2 <- ((Talla[j] + xs) - Lage[i]) / sd_age[i]
    pdf[i, j] <- pnorm(z2) - pnorm(z1)
  }
  pdf[i, ] <- pdf[i, ] / sum(pdf[i, ]+1e-10)
}

return(pdf)
}


# Función YPR-----------------------------------------------------------
per_recruit <- function(params, data_list) {
  with(as.list(data_list), {
    par_vals <- params
    names(par_vals) <- c("L50", "slope", "Fcr", "Lr", "a0", "cv")

    Sel    <- 1 / (1 + exp(-log(19) * (Lage - par_vals["L50"]) / par_vals["slope"]))
    Mad    <- 1 / (1 + exp(-log(19) * (Lage - parbiol[6]) / (parbiol[7]-parbiol[6])))
    Wage   <- exp(parbiol[4]) * Lage^parbiol[5]
    
    Fref<- seq(0,max(c(3*M,1.1*par_vals["Fcr"])),M/50) 
    N<- rep(1,length(Lage))
    YPR  <- rep(0,length(Fref))
    BPR  <- rep(0,length(Fref))
    dtm  <- parbiol[8]
    h    <- parbiol[9]
    target=parbiol[10]
    Fcur=par_vals["Fcr"]
    
    for (i in 2:Nages) {
      N[i]=N[i-1]*exp(-M)
    }
    N[i]=N[i]/(1-exp(-M))
    B0<- sum(N*Mad*Wage*exp(-dtm*M))
    alfa <- 4 * h / (5 * h - 1)
    beta <- (1 - h) / (5 * h - 1) * B0
    
    
    for (j in 1:length(Fref)) {
      F=Fref[j]*Sel
      Z=F+M
      
      for (i in 2:Nages) {
        
        N[i]=N[i-1]*exp(-Z[i-1])
      }
      N[i]=N[i]/(1-exp(-Z[i]))
      
      C<- N*F/Z*(1-exp(-Z))
      BPR[j]<- alfa*sum(N*Mad*Wage*exp(-dtm*Z))-beta
      YPR[j]<- sum(C*Wage)*alfa*BPR[j]/(beta+BPR[j])
      
      
      if(Fref[j]<Fcur){
        Ncurr=N
        Ccurr=C
        BPRcur=BPR[j]
        YPRcur=YPR[j]
      }
      
      if(j>1 && BPR[j]/BPR[1]>0.99*tar){
        Ftar=Fref[j]
        Ntar=N
        Ctar=C
        BPRtar=BPR[j]
        YPRtar=YPR[j]
      }
      
      if(j==1){
        N0=N
      }
    }
    
    
    apr_out=list(Fref=Fref,BPR=BPR,YPR=YPR,Ncurr=Ncurr,Ccurr=Ccurr,B0=B0,N0=N0,
                 Ntar=Ntar,Ctar=Ctar,Ftar=Ftar,B0=BPR[1],BPRtar=BPRtar,BPRcur=BPRcur,
                 YPRtar=YPRtar,YPRcur=YPRcur,Fcur=Fcur,Tmax=Tmax,edad=ages,Lage=Lage,
                 sd_age=sd_age)
    
    return(apr_out)
  })}


# Función Gráficos----------------------------------------
LBPA_graph=function(data_graph){
  with(as.list(data_graph), {

    nedades <-length(edad)

    pdf=pdf_mat(Talla,Lage, sd_age)

    Cagelength<-pdf
    Ctarlength<-pdf
    Nagelength<-pdf
    Nage0length<-pdf
    Ntarlength<-pdf

    for (i in 1:nedades) {
      Cagelength[i,]<-Ccurr[i]*pdf[i,]
      Ctarlength[i,]<-Ctar[i]*pdf[i,]
      Ntarlength[i,]<-Ntar[i]*pdf[i,]
      Nagelength[i,]<-Ncurr[i]*pdf[i,]
      Nage0length[i,]<-N0[i]*pdf[i,]
    }
    
    
    par(mfrow = c(1, 1))
    
    matplot(Talla,pobs,type="b",lty=1, col="darkgray",xlab="Length", ylab="Proportion",main="LBPA model fit",cex.main = 1.,lwd=2.)
    lines(Talla,ppred ,type="l", col="red",lwd=2)
    legend("topright",c("data","model"),col=c("gray","red"),
           lty=1,lwd=2,bty="n")
    abline(v=parbiol[1],col="blue",lty=4); text(parbiol[1],0.1,expression(L["00"]))
    abline(v=parbiol[6],col="blue",lty=4); text(parbiol[6],0.1,expression(L["50m"]))
    
    
    par(mfrow = c(1, 2))
    plot(Talla,ppred ,type="l", col="red",lwd=2,xlab="Length", ylab="Proportion",
         main="LBPA model fit", ylim=c(0,max(pobs)))
    
    
    for (i in 2:length(pobs[1,])) {
      lines(Talla,pobs[,i-1],type="p",col="gray",pch=20,cex=1.5)
    }
    
    lines(Talla,ppred, col="red",lwd=2)
    legend("topright",c("data","model"),col=c("gray","red"),
           lty=1,lwd=2,bty="n")
    
    
    unos=matrix(1,length(pobs[1,]),1)
    res=(pobs-ppred%*%t(unos))
    res=res/sd(res)
    
    hist(res,prob=T,main="Residuals",ylab="Density",xlab="Normalized residual")
    x <- seq(min(res), max(res), length = 200)
    lines(x, dnorm(x), col = "red", lwd = 2)
    box()

    
    x=seq(1,Tmax)
    Lx=parbiol[1]*(1-exp(-parbiol[2]*(x+.5)))
    
    par(mfcol = c(2, 2))
    plot(x,Lx,type="b",ylim=c(0,parbiol[1]), xlim=c(0,Tmax), xlab="Relative age",ylab="Length",
         main="Growth curve",lwd=2)
    lines(tr,Lage[1],pch=20,type="p",col="green",cex=3)
    grid()
    abline(h=parbiol[1],lty=2)
    abline(h=Lage[1],lty=2)
    abline(v=tr,lty=2)
    text(tr+1,Lage[1]*.9,paste("Lr= ",round(Lage[1],1)),col="red")
 
    L50=params[1];slope=params[2];L50m=parbiol[6];L95m=parbiol[7]
    
    plot(Talla,1/(1+exp(-log(19)*(Talla-L50)/(slope))),type="l",col="green",lwd=2,
         xlab="Length", ylab="Proportion",main="Selectivity and maturity",ylim=c(0,1))
    lines(Talla,1/(1+exp(-log(19)*(Talla-L50m)/(L95m-L50m))),type="l",col="blue",lwd=2)
    abline(h=0.5,lty=2)
    abline(v=L50,lty=2)
    legend("topright",c("Maturity","Selectivity"),col=c("blue","green"),
           lty=1,lwd=2,bty="n")
    text(L50*1.05,0.05,paste("L50=",round(L50,2)),col="red")
    
    
    plot(Fref,YPR_eq/max(YPR_eq),type="l", col="blue",lwd=2,xlab="Fishing mortality",ylab="Biomass, Yield",
         main="Equilibrium curves",ylim=c(0,1))
    
    if(Fcr>max(Fref)){text(max(Fref)*.95,0.05,paste("Fcr>",round(max(Fref),2)),col="red")}
    
    lines(Fref,BPR_eq/max(BPR_eq),type="l", col="green",lwd=2)
    abline(h=target,lty=2)
    abline(v=Ftar,lty=2)
    abline(v=Fcr,lty=2)
    
    lines(Fcr,SPR,type="p",pch=20)
    text(Ftar,target*1.1,paste("Fmsy=",round(Ftar,2)),col="red")
    text(Fcr*1.1,0.05,paste("Fcr=",round(Fcr,2)),col="red")
    legend("topright",c("Yield","Biomass"),col=c("blue","green"),
           lty=1,lwd=2,bty="n")
    
    
    eje=c("a) Unfished","b) Target","c) Current")
    resp=c(1,target,SPR)
    barplot(resp~eje,ylab="Proportion",xlab="Biomass",main="Proportion of B0",ylim=c(0,1.1),
            col=c("lightgreen","lightblue","gray"))
    abline(h=target,lty=2)
    box()
    text(3,SPR+0.05,paste("SPR=",round(SPR,2)),col="red")
    
    # 
    barplot(rbind(N0,Ncurr),beside=TRUE, width=1, space=c(0,0.2), names.arg = paste(edad),col=c("lightblue","pink"),xlab="Relative age",
            ylab="Density",border=NA,
            main="Population at-age")
    box()
    legend("topright",c("Unfished","Current"),col=c("lightblue","pink"),
           lty=1,lwd=3, bty="n")

    barplot(rbind(Ccurr,Ctar), beside=TRUE, width=1, space=c(0,0.2),border=NA,
            names.arg = paste(edad), col=c("lightblue","pink"),
            xlab="Relative age",ylab="Density",main="Catch at-age")
    legend("topright",c("Current","Target"),col=c("lightblue","pink","blue"),
           lty=1,lwd=2, bty="n")
    box()


    #
    areaplot(Talla,t(pdf)%*%Ctar,
             col="pink",
             xlim = c(min(Talla),1.1*max(Talla)),
             ylim = c(0,max(t(pdf)%*%Ctar,t(pdf)%*%Ccurr)),
             ylab="Density",
             xlab="Length",
             main="Catch at-length")
  
    lines(Talla,t(pdf)%*%Ccurr, col="black", lwd=2, lty=1)
    matlines(Talla,t(Cagelength),lty=3,col="black")
    Ltar=round(sum((t(pdf)%*%Ctar)*Talla)/sum(t(pdf)%*%Ctar),2)
    Lcur=round(sum((t(pdf)%*%Ccurr)*Talla)/sum(t(pdf)%*%Ccurr),2)

        

     legend("topright",c("Current","Target"),col=c("black","pink"),
            lty=1,lwd=2, bty="n")


    #kobe plot
    par(mfrow = c(1, 1))
    maxY=max(2,Fcr/Ftar*1.2)
    BPRtar=BPR_eq[1]*target
    maxX=max(2,BPR_eq/BPRtar)
    plot(0,0,type="l", col="gray",lwd=2,
         main="Kobe plot",
         ylab="F/Fmsy",xlab="B/Bmsy",ylim=c(0,maxY),xlim=c(0,maxX))
    
    polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1")
    polygon(c(1,maxX,maxX,1),c(0,0,1,1),col="green")
    polygon(c(1,maxX,maxX,1),c(1,1,maxY,maxY),col=rgb(1, 0.84, 0))
    polygon(c(0,1,1,0),c(1,1,maxY,maxY),col="tomato1")
    # 
    lines(SPR/target,Fcr/Ftar,type="p",pch=3,cex=3.0,lwd=2)
    lines(BPR_eq/BPRtar,Fref/Ftar,lty=2)
    text(SPR/target,Fcr/Ftar+0.1,paste("B/Bmsy=",round(SPR/target,2),"; F/Fmsy=",round(Fcr/Ftar,2)), font=2)
    
    
    
    
  })  
}

# Optimización de parámetros
pars_fin <- optim(par = parini, fn = ferrLBPA, data = data_list, method = "BFGS", hessian = TRUE)
mcov     <- solve(pars_fin$hessian)
parfin   <- exp(pars_fin$par)
sd_par   <- sqrt(diag(mcov)) * parfin


# Resultados
data_list$flag=2
r1  <-ferrLBPA(log(parfin), data_list)
r2 <- per_recruit(parfin, r1)


# Reporte de resultados resumido
table1 <- matrix(cbind(round(parfin, 2), round(sd_par, 3)), ncol = 2)
rownames(table1) <- c("Selectivity length at 50% (L50)", "Slope (d)","Fishing mortality (Fcr)",
                      "Size of recruits (Lr)","Invariant std in length (a0)", "Coeff of variation length at-age (cv)")
colnames(table1)<-c("Value","sd")

SPR=r2$BPRcur/r2$B0
table2 <- matrix(ncol=1, round(c(r2$B0, r2$BPRcur, r2$BPRtar, SPR, SPR/target, r2$Ftar,r2$Fcur/r2$Ftar,r2$YPRcur,r2$YPRtar),2))
rownames(table2) <- c("Virginal biomass per-recruit (BPR0)", "Current BPR", "Target BPR","Current spawning potential ratio (SPR)",
                      "Overexploitation index (SPR/SPRtar)", "Target fishing mortality (Ftar)","Overfishing index (F/Ftar)",
                      "Current yield per-recruit (YPRcur)","Target  yield per-recruit (YPRtar)")
colnames(table2)<-"Value"

table3 <- matrix(ncol=1, round(c(r1$lprop,r1$lprior,r1$lprop+sum(r1$lprior)),2))
rownames(table3) <- c("LF Proportions", "L50", "d","Fcr","Lr", "a0","cv","Total")
colnames(table3)<-"log-likelihood"


output <- list(table1 = table1, table2=table2, table3=table3, vars1=r1, vars2=r2)

data_graph=list(params=parfin,
                parbiol=parbiol,
                SPR=SPR,
                Talla=r1$talla,
                pobs=r1$pobs,
                ppred=r1$ppred,
                Ccurr=r2$Ccurr,
                Ctar=r2$Ctar,
                Ntar=r2$Ntar,
                Ncurr=r2$Ncurr,
                N0=r2$N0,
                Fref=r2$Fref,
                YPR_eq=r2$YPR,
                BPR_eq=r2$BPR,
                target=target,
                Ftar=r2$Ftar,
                Fcr=r2$Fcur,
                Tmax=r2$Tmax,
                edad=r2$edad,
                Lage=r2$Lage,
                sd_age=r2$sd_age,
                Lr=parfin[4],
                tr=r1$tr)
  
  if(graph_opt[1]==T){
    LBPA_graph(data_graph)
  }
  
  if(save_opt[1]==T){
    
    library(openxlsx)
    wb <- createWorkbook()
    addWorksheet(wb, "Table1_Parameters")
    addWorksheet(wb, "Table2_Variables")
    addWorksheet(wb, "Table3_Likelihood")
    
    
    writeData(wb, sheet = "Table1_Parameters", x = table1, rowNames = T )
    writeData(wb, sheet = "Table2_Variables", x = table2, rowNames = T )
    writeData(wb, sheet = "Table3_Likelihood", x = table3, rowNames = T )
    
    saveWorkbook(wb,paste("Outcomes_",name), overwrite = TRUE)
  } 
    

  return(output)
}

