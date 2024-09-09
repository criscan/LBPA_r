# LBPA_r

# Length-Based Pseudo-cohort Analysis (LBPA) (https://doi.org/10.1016/j.fishres.2020.105810)

Stock status for many medium- and small-scale fisheries is unknown due, for example, to a lack of catch data and the absence of scientific observer programs. However, length-frequency data are often available for such fisheries because they are the cheapest and easiest data to obtain. Several stock assessment methods have been developed that use length-frequency data and make equilibrium assumptions regarding both recruitment and fishing mortality. These assumptions raise questions regarding the reliability of the results, particularly when the method is applied to a single sample of length-frequency. We developed a Length-Based Pseudo-cohort Analysis (LBPA) model whose parameters can be estimated using multiple length frequencies and penalized maximum likelihood, under the assumption that using more than one length-frequency sample reduces the effects of the equilibrium conditions assumed in the model. This work provides guidelines that should be considered when using length-based pseudo-cohort models for data-poor fisheries.

The model is a package for R and can be installed as

1. If devtools package is already installed
   
   install_github("https://github.com/criscan/LBPA_r")

3. If devtools package is not installed

   devtools::install_github("https://github.com/criscan/LBPA_r")

When running the model,  be sure the Excel library is installed and the data file in the same path. The Excel data file has four sheets    

a. Lenght frequencies of catches by year/sample (LF)     
b. Biological parameters    

   Loo	k	= Growth parameters,
   M	= Natural mortality rate,
   log_aw =condition coefficient of length-weight relationship (in log scale),
   bw = potential coefficient of length-weight relationship,
   L50m,	L95m,	dtm = length of maturity at 50 and 95 percent, and year fraction when spawn occurs,	
   h= steepness of B&H stock-recruitment relationships, 	
   SPRtarget= management objective as B0 fraction.

c. Start values for model's parameters and its coefficient of variation. Each model parameter is assumed to follow a lognormal distribution with mean mu=log(a priori value) and standard deviation sd=cv.
d. Weighting factors for each LF considered    

The model can be run as model1=LBPA_fits("filename.xlsx",graph_opt=T or F, save_opt=T or F). An example of a data file (.xlsx) was added and can be downloaded before package installation. After running, if save_opt=T, an Excel (xlsx) file will be produced with five sheets: model parameters, model variables, likelihood components values, per recruits analysis, and different population variables by age. In the same context, several graphics of the model will be generated if graph_opt=T

For any question: cristian.canales.r@pucv.cl
