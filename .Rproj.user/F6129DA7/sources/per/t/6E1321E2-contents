#' @title Length-Based Pseudo-cohort Analysis - LBPA
#' @description The model estimates the fishing mortality and depletion of the Spawning Potential Ratio (SPR) in an equilibrium condition.
#' The model is fitted to several length frequencies of catches simultaneously by minimizing a penalized log-likelihood function.
#' @param name excel (.xlsx) file with data and parameters. This file contains four sheets: LF data, biological parameters,
#' initial parameters and data weighting
#' @return Several graphics with de model fit and population variables are displayed. Four csv files are generated: model parameters,
#' population variables, log-likelihood and per-recruit variables
#' @examples
#' LBPA_fits("data_file.xlsx")


LBPA_fits=function(name){

#' @description The model estimates the fishing mortality and depletion of the Spawning Potential Ratio (SPR) in an equilibrium condition.
#' The model is fitted to several length frequencies of catches simultaneously by minimizing a penalized log-likelihood function.
#' @param name excel (.xlsx) file with data and parameters. This file contains four sheets: LF data, biological parameters,
#' initial parameters and data weighting
#' @return Several graphics with de model fit and population variables are displayed. Four csv files are generated: model parameters,
#' population variables, log-likelihood and per-recruit variables
#' @examples
#' LBPA_fits("data_file.xlsx")



  ferrLBPA<-function(parini,data_list){


    # Por recluta
    apr_out <- function(Sel, M, Mad, Wage, Tmax) {

      Fref=seq(0,3*M,M/50)
      N=rep(1,Tmax)
      YPR=rep(0,length(Fref))
      BPR=rep(0,length(Fref))

      for (j in 1:length(Fref)) {
        F=Fref[j]*Sel
        Z=F+M

        for (i in 2:Tmax) {

          N[i]=N[i-1]*exp(-Z[i-1])
        }
        N[i]=N[i]/(1-exp(-Z[i]))

        C=N*F/Z*(1-exp(-Z))
        YPR[j]=sum(C*Wage)
        B0=sum(N0*Mad*Wage*exp(-dtm*M))
        BPR[j]=sum(N*Mad*Wage*exp(-dtm*Z))

      }
      apr_out=list(Fref=Fref,BPR=BPR,YPR=YPR)

      return(apr_out)
    }



    Loo=data_list$parbiol[1]
    k=data_list$parbiol[2]
    M=data_list$parbiol[3]
    ln_a=data_list$parbiol[4]
    b=data_list$parbiol[5]
    L50m=data_list$parbiol[6]
    L95m=data_list$parbiol[7]
    dtm=data_list$parbiol[8]
    h=data_list$parbiol[9]
    target=data_list$parbiol[10]
    nm=data_list$nm
    lambda=data_list$lambda

    Tmax=round(-log(0.025)/M)
    age=rep(0,Tmax)
    Lage=rep(0,Tmax)

    datos=data_list$datos
    Talla=datos[,1]
    fobs=datos[,2:length(datos[1,])]

    L50=exp(parini[1])
    slope=exp(parini[2])
    Fcr=exp(parini[3])
    Lr=exp(parini[4])
    a0=exp(parini[5])
    cv=exp(parini[6])


    N=rep(1,Tmax)
    N0=rep(1,Tmax)
    Ntar=rep(1,Tmax)
    Ctar=rep(1,Tmax)



    Lage[1]=Lr
    for (i in 2:Tmax) {
      Lage[i]=Lage[i-1]*exp(-k)+Loo*(1-exp(-k))
    }

    sd_age=a0+cv*Lage

    Sel=1/(1+exp(-log(19)*(Lage-L50)/(slope)))
    Mad=1/(1+exp(-log(19)*(Lage-L50m)/(L95m)))
    Wage=exp(ln_a)*Lage^b

    Sel_l=1/(1+exp(-log(19)*(Talla-L50)/(slope)))


    F=Fcr*Sel
    Z=F+M

    for (i in 2:Tmax) {

      N0[i]=N0[i-1]*exp(-M)
      N[i]=N[i-1]*exp(-Z[i-1])
    }
    N0[i]=N0[i]/(1-exp(-M))
    N[i]=N[i]/(1-exp(-Z[i]))

    C=N*F/Z*(1-exp(-Z))
    B0=sum(N0*Mad*Wage*exp(-dtm*M))
    B=sum(N*Mad*Wage*exp(-dtm*Z))

    R0=1
    alfa=4*h*R0/(5*h-1);
    beta=(1-h)/(5*h-1)*B0;

    SPR=(alfa*B-beta)/B0

    # Matriz clave talla-edad inversa

    dtalla=Talla[2]-Talla[1]
    pdf=matrix(NA,Tmax,length(Talla))
    xs=0.5*dtalla

    for(i in 1: Tmax) #loop over ages
    {
      for(j in 1: length(Talla)) #loop over ages
      {
        z1=((Talla[j]-xs)-Lage[i])/sd_age[i]
        z2=((Talla[j]+xs)-Lage[i])/sd_age[i]
        pdf[i,j]=pnorm(z2)-pnorm(z1)}
      pdf[i,]=pdf[i,]/sum(pdf[i,])

    }

    Ctalla=t(pdf)%*%C


    ppred=Ctalla/sum(Ctalla)
    sum=0

    for (i in 2:length(datos[1,])) {
      sum=sum+nm*lambda[i-1]*sum(fobs[,i-1]/sum(fobs[,i-1])*log(ppred))
    }

    #penalizacion de parametros
    penal=0.5*sum((data_list$prioris-log(c(L50,slope,Fcr,Lr,a0,cv)))^2/data_list$cv_par^2)

    fun=-sum+penal


    return(fun)


  }


  grafLBPA<-function(parfin,data_list){


    # Por recluta
    apr_out <- function(Sel, M, Mad, Wage, Tmax, Fcr) {

      Fref=seq(0,3*M,M/50)
      N=rep(1,Tmax)
      YPR=rep(0,length(Fref))
      BPR=rep(0,length(Fref))

      for (j in 1:length(Fref)) {
        F=Fref[j]*Sel
        Z=F+M

        for (i in 2:Tmax) {

          N[i]=N[i-1]*exp(-Z[i-1])
        }
        N[i]=N[i]/(1-exp(-Z[i]))

        C=N*F/Z*(1-exp(-Z))
        YPR[j]=sum(C*Wage)
        BPR[j]=sum(N*Mad*Wage*exp(-dtm*Z))

        if(Fref[j]<Fcr){
          Ncurr=N
          Ccurr=C
        }

        if(j==1){
          N0=N
        }

      }
      B0=sum(N0*Mad*Wage*exp(-dtm*M))
      apr_out=list(Fref=Fref,BPR=BPR,YPR=YPR,Ncurr=Ncurr,Ccurr=Ccurr,B0=B0,N0=N0)

      return(apr_out)
    }



    Loo=data_list$parbiol[1]
    k=data_list$parbiol[2]
    M=data_list$parbiol[3]
    ln_a=data_list$parbiol[4]
    b=data_list$parbiol[5]
    L50m=data_list$parbiol[6]
    L95m=data_list$parbiol[7]
    dtm=data_list$parbiol[8]
    h=data_list$parbiol[9]
    target=data_list$parbiol[10]
    nm=data_list$nm
    lambda=data_list$lambda

    Tmax=round(-log(0.025)/M)
    age=seq(1,Tmax)
    Lage=rep(0,Tmax)


    datos=data_list$datos
    Talla=datos[,1]
    fobs=datos[,2:length(datos[1,])]

    L50=parfin[1]
    slope=parfin[2]
    Fcr=parfin[3]
    Lr=parfin[4]
    a0=parfin[5]
    cv=parfin[6]


    N=rep(1,Tmax)
    N0=rep(1,Tmax)
    Ntar=rep(1,Tmax)
    Ctar=rep(1,Tmax)



    Lage[1]=Lr
    for (i in 2:Tmax) {
      Lage[i]=Lage[i-1]*exp(-k)+Loo*(1-exp(-k))
    }

    sd_age=a0+cv*Lage

    Sel=1/(1+exp(-log(19)*(Lage-L50)/(slope)))
    Mad=1/(1+exp(-log(19)*(Lage-L50m)/(L95m)))
    Wage=exp(ln_a)*Lage^b

    Sel_l=1/(1+exp(-log(19)*(Talla-L50)/(slope)))

    ypr <- apr_out(Sel, M, Mad, Wage, Tmax, Fcr)

    B0=ypr$B0
    N0=ypr$N0
    C=ypr$Ccurr
    N=ypr$Ncurr

    Z=Fcr*Sel+M
    B=sum(N*Mad*Wage*exp(-dtm*Z))


    R0=1
    alfa=4*h*R0/(5*h-1);
    beta=(1-h)/(5*h-1)*B0;

    SPR=(alfa*B-beta)/B0

    # Matriz clave talla-edad inversa

    dtalla=Talla[2]-Talla[1]
    pdf=matrix(NA,Tmax,length(Talla))
    xs=0.5*dtalla

    for(i in 1: Tmax) #loop over ages
    {
      for(j in 1: length(Talla)) #loop over ages
      {
        z1=((Talla[j]-xs)-Lage[i])/sd_age[i]
        z2=((Talla[j]+xs)-Lage[i])/sd_age[i]
        pdf[i,j]=pnorm(z2)-pnorm(z1)}
      pdf[i,]=pdf[i,]/sum(pdf[i,])

    }


    Ctalla=t(pdf)%*%C


    ppred=Ctalla/sum(Ctalla)
    sum=0

    for (i in 2:length(datos[1,])) {
      sum=sum+nm*lambda[i-1]*sum(fobs[,i-1]/sum(fobs[,i-1])*log(ppred))
    }

    #penalizacion de parametros
    penal=sum(0.5*(data_list$prioris-log(c(L50,slope,Fcr,Lr,a0,cv)))^2/data_list$cv_par^2)


    par(mfrow = c(2, 2))

    pobs=fobs/matrix(1,length(datos[,1]),1)%*%colSums(fobs)

    matplot(Talla,pobs,type="l",lty=1, col="gray",xlab="Talla", ylab="Frecuencia",main="Frec. tallas",cex.main = 1.)


    plot(Talla,ppred ,type="l", col="red",lwd=2,xlab="Talla", ylab="Frecuencia",
         main="Frec. tallas", ylim=c(0,max(pobs)),cex.main = 1.)
    pobs=fobs


    for (i in 2:length(datos[1,])) {
      lines(Talla,fobs[,i-1]/sum(fobs[,i-1]),type="p",col="gray",pch=20,cex=1.5)
      pobs[,i-1]=fobs[,i-1]/sum(fobs[,i-1])
    }
    lines(Talla,ppred, col="red",lwd=2)
    legend("topright",c("dato","modelo"),col=c("gray","red"),
           lty=1,lwd=2,bty="n",cex=0.8)



    unos=matrix(1,length(fobs[1,]),1)
    res=(pobs-ppred%*%t(unos))
    res=res/sd(res)

    hist(res,prob=T,main="Residuales",cex.lab = 1.,
         cex.axis = 1.,ylab="Densidad",xlab="Residual normalizado",cex.main = 1.)
    x <- seq(min(res), max(res), length = 200)
    lines(x, dnorm(x), col = "red", lwd = 2)
    box()

 par(mfrow = c(2, 2))

 pobs=fobs/matrix(1,length(datos[,1]),1)%*%colSums(fobs)

 matplot(Talla,pobs,type="l",lty=1, col="gray",xlab="Length", ylab="Frequency",main="Length frequency",cex.main = 1.)


 plot(Talla,ppred ,type="l", col="red",lwd=2,xlab="Length", ylab="Frequency",
      main="Length frequency", ylim=c(0,max(pobs)),cex.main = 1.)
 pobs=fobs


 for (i in 2:length(datos[1,])) {
  lines(Talla,fobs[,i-1]/sum(fobs[,i-1]),type="p",col="gray",pch=20,cex=1.5)
  pobs[,i-1]=fobs[,i-1]/sum(fobs[,i-1])
 }
 lines(Talla,ppred, col="red",lwd=2)
 legend("topright",c("data","model"),col=c("gray","red"),
        lty=1,lwd=2,bty="n",cex=0.8)



 unos=matrix(1,length(fobs[1,]),1)
 res=(pobs-ppred%*%t(unos))
 res=res/sd(res)

 hist(res,prob=T,main="Residuals",cex.lab = 1.,
      cex.axis = 1.,ylab="Density",xlab="Normalized residual",cex.main = 1.)
 x <- seq(min(res), max(res), length = 200)
 lines(x, dnorm(x), col = "red", lwd = 2)
 box()

 par(mfrow = c(2, 2))

 plot(Talla,pdf[1,],type="l",col="green",lwd=2,
      xlab="Length", ylab="Proportion",main="Recruitment and modal lengths",cex.main = 1.)
 abline(v=Lage,lty=2,col="gray")
 text(Lr*1.1,0.01,paste("Lr=",round(Lr,2)))


 Sel=1/(1+exp(-log(19)*(Lage-L50)/(slope)))
 Mad=1/(1+exp(-log(19)*(Lage-L50m)/(L95m)))

  plot(Talla,1/(1+exp(-log(19)*(Talla-L50)/(slope))),type="l",col="green",lwd=2,
       xlab="Length", ylab="Proportion",main="Selectivity and maturity",cex.main = 1.)
  lines(Talla,1/(1+exp(-log(19)*(Talla-L50m)/(L95m-L50m))),type="l",col="blue",lwd=2)
  abline(h=0.5,lty=2)
  abline(v=L50,lty=2)
  legend(max(Talla)*0.7,0.95,c("Mat","Select"),col=c("blue","green"),
         lty=1,lwd=2,bty="n",cex=0.8)
  text(L50*1.1,0.05,paste("L50=",round(L50,2)))


  plot(age,Lage,type="b",ylim=c(0,Loo),xlab="Relative age", ylab="Length", main="Growth model age-length", cex.main=1,
       lwd=2)
  abline(h=Lage,lty=2,col="gray")


  par(mfrow = c(2, 1))

  ypr <- apr_out(Sel, M, Mad, Wage, Tmax, Fcr)

  BPR_eq=alfa*ypr$BPR-beta
  Reclu=alfa*BPR_eq/(beta+BPR_eq)
  YPR_eq=ypr$YPR*Reclu


  id=which(BPR_eq/BPR_eq[1]<target)[1]

  plot(ypr$Fref,YPR_eq/max(YPR_eq),type="l", col="blue",lwd=2,xlab="Fishing mortality",ylab="Biomass, Yield",
       main="Equilibrium curves",ylim=c(0,1),cex.main = 1.)
  lines(ypr$Fref,BPR_eq/max(BPR_eq),type="l", col="green",lwd=2)
  abline(h=target,lty=2)
  abline(v=ypr$Fref[id],lty=2)
  abline(v=Fcr,lty=2)

  lines(Fcr,SPR,type="p",pch=20,cex=2.0)
  text(ypr$Fref[id]*1.2,0.1,paste("Fmsy=",round(ypr$Fref[id],2)),col="red",cex=0.8)
  text(Fcr*1.1,0.1,paste("Fcr=",round(Fcr,2)),col="red",cex=0.8)
  legend("topright",c("Yield","Biomass"),col=c("blue","green"),
         lty=1,lwd=2,bty="n",cex=0.8)



  eje=c("a) Unfished","b) Target","c) Current")
  resp=c(1,target,SPR)
  barplot(resp~eje,ylab="Proportion",xlab="Biomass",main="Proportion of B0",ylim=c(0,1.1),
          col=c("lightgreen","lightblue","gray"),cex.main = 1.)
  abline(h=target,lty=2)
  box()
  text(3,1.5*SPR,paste("SPR=",round(SPR,2)),col="red",cex=0.8)



  YPRtar=YPR_eq[id]
  Ftar=ypr$Fref[id]
  id2=which(BPR_eq/BPR_eq[1]<SPR)[1]
  YPRcur=YPR_eq[id2]

  ypr <- apr_out(Sel, M, Mad, Wage, Tmax, Ftar)
  Ctar=ypr$Ccurr
  Ntar=ypr$Ncurr


  Cagelength<-pdf
  Ctarlength<-pdf
  Nagelength<-pdf
  Nage0length<-pdf
  Ntarlength<-pdf


  for (i in 1:Tmax) {
    Cagelength[i,]<-C[i]*pdf[i,]#/sum(C)
    Ctarlength[i,]<-Ctar[i]*pdf[i,]#/sum(C)
    Ntarlength[i,]<-Ntar[i]*pdf[i,]#/sum(C)
    Nagelength[i,]<-N[i]*pdf[i,]#/sum(C)
    Nage0length[i,]<-N0[i]*pdf[i,]#/sum(C)
  }

  par(mfrow = c(2,2))


  edad=c(1:Tmax)
  barplot(N0~edad,col="lightblue",xlab="Relative age",ylab="Density",
          main="Population at-age",cex.main = 1.)
  barplot(N~edad, add = T,col = "gray")
  barplot(C~edad, add = T,col = "blue")
  box()
  legend("topright",c("Unfished","Current","Catch"),col=c("lightblue","gray","blue"),
         lty=1,lwd=2, bty="n",cex=0.8)


    plot(Talla,colSums(Nagelength), type="l", lwd=2, col="blue",
       ylab="Density",
       xlab="Length",
       main="Population at-length",cex.main = 1.)

  lines(Talla,colSums(Ntarlength),
        type="l", cex.lab = 1.5,
        lwd=2, lty=2,
        col="black",
        xlim = c(min(Talla),1.1*max(Talla)),
        ylim = c(0,max(Ctarlength)))

  lines(Talla,colSums(Nage0length),
        type="l", cex.lab = 1.5,
        lwd=2, lty=2,
        col="green",
        xlim = c(min(Talla),1.1*max(Talla)),
        ylim = c(0,max(Nagelength)))



  legend("topright",c("Current","Target","Unfished"),col=c("blue","black","green"),
         lty=c(1,2,1),lwd=2,bty="n",cex=0.8)


  lines(Talla,Nagelength[1,],
        type="l", cex.lab = 1.5,
        xlim = c(min(Talla),1.1*max(Talla)),
        ylim = c(0,max(Nagelength)))



  for (i in 2:Tmax) {
    lines(Talla,Nagelength[i,], type="l")
  }


  barplot(C~edad,col="lightblue",xlab="Relative age",ylab="Density",
          main="Catch at-age",cex.main = 1.)
  box()



  plot(Talla,ppred, type="l", lwd=2, col="blue",
       ylab="Density",
       xlab="Length",
       main="Catch at-length",cex.main = 1.)

  Ctarlength=Ctarlength/sum(Ctarlength)
  lines(Talla,colSums(Ctarlength),
        type="l", cex.lab = 1.5,
        lwd=2, lty=2,
        col="black",
        xlim = c(min(Talla),1.1*max(Talla)),
        ylim = c(0,max(Ctarlength)))

  Cagelength=Cagelength/sum(Cagelength)
  lines(Talla,Cagelength[1,],
        type="l", cex.lab = 1.5,
        xlim = c(min(Talla),1.1*max(Talla)),
        ylim = c(0,max(Cagelength)))


  for (i in 2:Tmax) {
    lines(Talla,Cagelength[i,], type="l")
  }

  legend("topright",c("Current","Target"),col=c("blue","black"),
         lty=c(1,2),lwd=2,bty="n",cex=0.8)





  like=c(-sum,0.5*(data_list$prioris-log(c(L50,slope,Fcr,Lr,a0,cv)))^2/data_list$cv_par^2)


  outputs=list(edad=age,talla=Lage,sd_edad=sd_age,Select=Sel,Mad=Mad,Peso=Wage,N0=N0,N=N,F=F,Z=Z,C=C,alfa=alfa,beta=beta,ptalla=pdf,
               Fref=ypr$Fref,BPReq=BPR_eq,YPReq=YPR_eq,pR0=Reclu,SPR=SPR,Fcur=Fcr,Ftar=Ftar,YPRtar=YPRtar,YPRcur=YPRcur,likeval=like)

    return(outputs)


  }



  graphics.off()
  library(readxl)

  # Archivo con los datos

  datos=as.matrix(read_xlsx(name,sheet=1,col_names = TRUE))
  parbiol=as.matrix(read_xlsx(name,sheet=2,col_names = TRUE))
  pars_ini=as.matrix(read_xlsx(name,sheet=3,col_names = TRUE))
  lambda=as.matrix(read_xlsx(name,sheet=4,col_names = TRUE))

  cvpar=pars_ini[2,]
  Talla=datos[,1]
  M=parbiol[3]


  # valores iniciales
  L50= pars_ini[1,1]
  slope=pars_ini[1,2]
  Fcr=pars_ini[1,3]
  Lr=pars_ini[1,4]
  a0=pars_ini[1,5]
  cv=pars_ini[1,6]


  parini=log(c(L50,slope,Fcr,Lr,a0+1e-5,cv+1e-5))


#tamaño muestra
nm=250

data_list=list(datos=datos,parbiol=parbiol,nm=nm,cv_par=cvpar,prioris=parini,lambda=lambda)
pars_fin=optim(par=parini,fn=ferrLBPA, data=data_list, method="BFGS")
target=parbiol[10]


parfin=exp(pars_fin$par)
solucion=data.frame(L50=parfin[1],slope=parfin[2],Fcr=parfin[3],Lr=parfin[4],a0=parfin[5],cv=parfin[6],LL=pars_fin$value)


v=grafLBPA(parfin,data_list)


table1 <- matrix(ncol=1, round(solucion[1:6], 2))
rownames(table1) <- c("Selectivity length at 50% (L50)", "Slope (d)","Fishing mortality (Fcr)",
                      "Size of recruits (Lr)","Invariant std in length (a0)", "Coeff of variation length at-age (cv)")

B0=v$BPReq[1]

 table2 <- matrix(ncol=1, round(c(B0,B0*v$SPR, B0*target, v$SPR ,target, v$Ftar,v$Fcur/v$Ftar,v$YPRcur,v$YPRtar),2))
 rownames(table2) <- c("Virginal biomass per-recruit (BPR0)", "Current BPR", "Target BPR","Current spawning potential ratio (SPR)",
                         "Target SPR (SPRtar)", "Target fishing mortality (Ftar)","Overfishing index (F/Ftar)",
                         "Current yield per-recruit (YPRcur)","Target  yield per-recruit (YPRtar)")

 table3 <- matrix(ncol=1, round(c(v$likeval,pars_fin$value),2))
 rownames(table3) <- c("LF Proportions", "L50", "d","Fcr","Lr", "a0","cv","Total")

 table4=data.frame(Fref=v$Fref,BPReq=v$BPReq,YPReq=v$YPReq)

 write.csv(table1, 'Parameters_LBPA.csv',row.names = T)

 write.csv(table2, 'Variables_LBPA.csv', row.names = T)

 write.csv(table3, 'logLL_LBPA.csv', row.names = T)

 write.csv(table4, 'Per_recruit.csv', row.names = F)



 print(knitr::kable(table1,"simple",caption = "1: Estimated parameters"))
 print(knitr::kable(table2,"simple",caption = "2: Per-recruit population's variables"))
 print(knitr::kable(table3,"simple",caption = "3: Log-likelihood components"))


}
