mPLPrd R CODE
Version 1.04.1, 21 May 2026

Authors:
Giuseppe Cavaliere, Silvia Goncalves, Morten O. Nielsen, and Edoardo Zanelli


1. Overview
-----------

mPLPrd is research code for regression-discontinuity inference using the
modified prepivoted local-polynomial procedure described in the associated
CGNZ paper, "Improved inference for nonparametric regression and
regression-discontinuity designs".

This folder contains the R implementation for version 1.04.1. It is a
script-style R package: source mPLP_functions.R first and mPLPrd.R second.

Version 1.04 is based on version 1.03 and adds an optional fixed-regressor
wild bootstrap confidence interval. Version 1.04.1 also restores the original
mplprdrobust() mPLP bias-correction window: auxiliary local fits for the mPLP
bias correction are computed on the 2h neighborhood, matching the paper's
mPLP/PLP support calculations and the simulation code used for the paper.


2. Files
--------

mPLPrd.R
  Main user-facing function. It computes the RD estimate, mPLP bias-corrected
  estimate, standard error, normal confidence interval, and optionally the
  fixed-regressor wild bootstrap confidence interval. The mPLP bias correction
  uses auxiliary local fits on |x_i - c| < 2h, as in the original
  mplprdrobust() simulation implementation.

mPLP_functions.R
  Auxiliary functions for kernels, local-polynomial estimation, equivalent
  kernels, mPLP standard errors, CCT-style standard errors, and the simplified
  local-polynomial standard error used by the bootstrap prepivoting step.

CGNZ_SRD_parall_prova.R
  Original-style example/simulation script copied forward from earlier
  versions.

CGNZ_SRD_bootstrap_rdrobust_bwselect_mc_1.04.1.R
  Monte Carlo script comparing rdrobust, the analytic mPLP interval, and the
  optional bootstrap interval using rdrobust bandwidth selectors.

CHANGES_1.04.txt
  Short summary of changes from version 1.03 to version 1.04.


3. Quick Start
--------------

  source("mPLP_functions.R")
  source("mPLPrd.R")

  fit <- mPLPrd(
    y = y,
    x = x,
    h = h,
    c = 0,
    p = 1,
    kersel = "tri",
    res = "cct-hc0",
    alpha = 0.05,
    fast = TRUE,
    g.loo = FALSE
  )

  fit$tau.hat
  fit$tau.hat_mlpbc
  fit$se_mplp
  fit$ci_mplp

Bootstrap confidence interval:

  fit.boot <- mPLPrd(
    y = y,
    x = x,
    h = h,
    c = 0,
    p = 1,
    kersel = "tri",
    res = "cct-hc0",
    alpha = 0.05,
    fast = TRUE,
    g.loo = FALSE,
    boot.ci = TRUE,
    boot.B = 999,
    boot.seed = 170498,
    boot.wild = "rademacher"
  )

  fit.boot$ci_mplp
  fit.boot$ci_mplp_boot


4. Main Function
----------------

  mPLPrd(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
         alpha = NULL, fast = NULL, g.loo = NULL,
         boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
         boot.wild = NULL, boot.keep = NULL)

New bootstrap arguments in version 1.04:

  boot.ci
    If TRUE, compute the fixed-regressor wild bootstrap confidence interval.
    Default: FALSE.

  boot.B
    Number of bootstrap repetitions. Default: 999 when boot.ci = TRUE.

  boot.seed
    Optional numeric seed set before drawing bootstrap weights.

  boot.wild
    Wild bootstrap weights. Accepted values are:
      "normal"      standard normal weights
      "rademacher"  two-point weights, -1 or 1
    Default: "rademacher".

  boot.keep
    If TRUE, return the bootstrap statistic draws. Default: TRUE.


5. Bootstrap Implementation
---------------------------

The bootstrap option implements the fixed-regressor wild bootstrap DGP in
equation (3.1) of the paper, adapted to the sharp RD estimator.

For observations inside the effective bandwidth window used by the cutoff
local-polynomial fits, version 1.04 computes side-specific local-polynomial
fitted values. If the internal bandwidth regularization enlarges the cutoff
window on one side, the bootstrap uses that enlarged side-specific window:

  ghat_+(x_i) for x_i >= c,
  ghat_-(x_i) for x_i < c.

Bootstrap outcomes are generated only for those observations:

  y_i^* = ghat_+(x_i) + epshat_i e_i,  x_i >= c and |x_i - c| <= h_+,
  y_i^* = ghat_-(x_i) + epshat_i e_i,  x_i <  c and |x_i - c| <= h_-,

where the regressors are fixed, epshat_i are the residuals already constructed
by the selected residual option, e_i are wild bootstrap weights, and h_+ and
h_- denote the effective right and left cutoff bandwidths after regularization.

For each bootstrap sample, the package recomputes the right and left
local-polynomial estimates at the cutoff and forms the RD-adapted modified
bootstrap statistic

  T_rd^* = Q.plus  * (ghat_+^*(c) - ghat_+(c))
         - Q.minus * (ghat_-^*(c) - ghat_-(c)).

This is the RD analogue of the bootstrap statistic in Section 3.2, with the
side-specific boundary scaling described for RD in Section 4.3.

The bootstrap confidence interval follows equation (2.6). The empirical
bootstrap cdf of T_rd^* is evaluated at the prepivoted probability levels
based on the inverse of equation (2.8):

  Hhat^{-1}(u) = Phi(mhat * Phi^{-1}(u)),
  mhat    = se_mplp / se_us,

where se_mplp is the mPLP standard error and se_us is the analogous
local-polynomial standard error for T_rd^*, including Q.plus and Q.minus,
but without the mPLP bias-correction component.

When boot.wild = "normal", T_rd^* is conditionally Gaussian. The package
still computes the interval by actual simulation of the fixed-regressor wild
bootstrap. Therefore ci_mplp_boot should approach the non-bootstrap mPLP
interval ci_mplp as boot.B grows, up to Monte Carlo error from the simulated
bootstrap quantiles.


6. Returned Bootstrap Values
----------------------------

When boot.ci = TRUE, the returned list also includes:

  ci_mplp_boot
    Fixed-regressor wild bootstrap confidence interval.

  boot.B
    Requested number of bootstrap repetitions.

  boot.success
    Number of finite bootstrap statistic draws.

  boot.probs
    Prepivoted probability levels used in the bootstrap quantiles.

  boot.mhat
    mhat = se_mplp / se_us from equation (2.8).

  boot.wild
    Wild-bootstrap weight distribution used.

  se_us
    Local-polynomial standard error for T_rd^*, including Q.plus and Q.minus,
    but without the mPLP bias-correction component.

  boot.t
    Bootstrap statistic draws, returned only when boot.keep = TRUE.


7. Verification
---------------

Latest checks on this machine:

  Patched mPLPrd() matches the old mplprdrobust() output up to numerical
  precision on identical simulated RD datasets.

Run:

  MPLPRD_MC_REP=2500 MPLPRD_BOOT_B=2 MPLPRD_MC_FAST=FALSE \
  MPLPRD_MC_G_LOO=FALSE MPLPRD_ERROR_DGP=gaussian MPLPRD_BW_VCE=nn \
  MPLPRD_BWSELECT_GRID=cerrd \
  Rscript CGNZ_SRD_bootstrap_rdrobust_bwselect_mc_1.04.1.R

Gaussian DGP1, n = 500, CCF-CE bandwidth:

  average bandwidth: 0.07334182
  rdrobust coverage: 0.936
  analytic mPLPrd coverage: 0.952
  analytic mPLPrd average length: 0.5416621
  rdrobust average length: 0.6801869


8. Notes
--------

- The bootstrap option is computationally more expensive than the analytic
  mPLP interval because it recomputes the cutoff estimator for every bootstrap
  repetition.

- The kernel aliases from version 1.03 are retained: "tri"/"triangular",
  "epa"/"epanechnikov", and "uni"/"uniform" work consistently.

- The inherited bwcheck = 21 bandwidth regularization remains.

- This is research replication code, not a formal CRAN-style package with a
  DESCRIPTION and NAMESPACE.
