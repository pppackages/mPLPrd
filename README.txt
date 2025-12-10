######################################################################################
## mPLPrd: an R package for efficient RDD inference
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
######################################################################################

mPLPrd performs estimation and inference in sharp regression-discontinuity designs by 
means of the mPLP method described in Cavaliere et al. (2025). 

This method estimates the ATE via local polynomial estimation at both sides of the cutoff,
and bias corrects by means of a local polynomial bootstrap. By prepivoting this bootstrap
DGP, robust bias correct standard errors are obtained.

