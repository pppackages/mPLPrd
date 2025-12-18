######################################################################################
## mPLPrd: an R package for efficient RDD inference
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
######################################################################################

mPLPrd performs estimation and inference in sharp regression-discontinuity designs by 
means of the mPLP method described in Cavaliere et al. (2025). 

This method estimates the ATE via local polynomial estimation at both sides of the cutoff,
and bias corrects by means of a local polynomial bootstrap. By prepivoting this bootstrap
DGP, robust bias correct standard errors are obtained.

########## Usage ###########

mPLPrd <- function(y, x, h, c = NULL, p = NULL, kersel = NULL, 
                  res = NULL, alpha = NULL, fast = NULL, g.loo = NULL) 


########## Arguments ###########


y            is the dependent variable 

x            is the running variable (or forcing variable, or score)

h            is the bandwidth parameter (to be chosen externally)

c            is the cutoff value. Default is c=0.

p            is the polynomial order used for estimating the local polynomials at 
             the left and at the right of the cutoff. It can take any odd and positive
             value. Default is p=1.

kersel       selects the kernel function used for the local polynomial estimators.
             The possible options are "epa" (epanechnikov's), "uni" (uniform) and 
             "tri" (triangular) kernel. Default is kersel="tri".

res          selects the residual class used to evaluate the standard erros. 
             Options are "loo", "cct-hc0", "cct-hc1", "cct-hc2", "cct-hc3", where the first
             one uses as fitted value the nonparametric regression curve at the right and 
             at the left of the cutoff (estimated at each value of x and using the polynomial 
             order p); while the remaining four options use the heteroskedasticity-robust residuals
             in Calonico et al. (2014, Econometrica), with weights 0, 1, 2 and 3, respectively.

alpha        is the significance level used for the confidence intervals.




