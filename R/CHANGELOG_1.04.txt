Changes of 1.04 mPLPrd from previous versions
==========================================


- The mPLP bias-correction auxiliary fits now use the 2h neighborhood,
  matching the original mplprdrobust() simulation implementation and the
  support of the mPLP/PLP bias-corrected weights in the paper.

- This corrects an unintended drift in the packaged mPLPrd() code where the
  auxiliary fits were computed over the h neighborhood unless res = "loo".
  That narrower window changed the mPLP bias correction and could materially
  distort finite-sample coverage.

- A direct check on identical simulated RD datasets confirms that patched
  mPLPrd() matches the old mplprdrobust() output up to numerical precision.

Version 1.04 adds an optional fixed-regressor wild bootstrap confidence
interval. The bootstrap option is not the default. Calls with boot.ci = FALSE
use the analytic mPLP interval; in 1.04.1 that interval follows the original
mplprdrobust() 2h auxiliary-window implementation.


New arguments
-------------

- boot.ci
  Set to TRUE to compute the bootstrap confidence interval.

- boot.B
  Number of bootstrap repetitions. Default: 999 when boot.ci = TRUE.

- boot.seed
  Optional seed for bootstrap weights.

- boot.wild
  Wild bootstrap weights: "normal" or "rademacher". Default: "rademacher".

- boot.keep
  If TRUE, return bootstrap statistic draws. Default: TRUE.

Bootstrap method
----------------

- Implements the fixed-regressor wild bootstrap DGP from equation (3.1) of Cavaliere, Goncalves, Nielsen & Zanelli (2026), adapted to the sharp RD setting.

- Bootstrap fitted values ghat(x_i) are computed only for observations inside
  the effective cutoff bandwidth window. If bandwidth regularization enlarges
  the window separately on the right or left, the bootstrap uses those
  side-specific effective windows.

- The bootstrap statistic is the RD-adapted modified statistic:

    T_rd^* = Q.plus  * (ghat_+^*(c) - ghat_+(c))
           - Q.minus * (ghat_-^*(c) - ghat_-(c)).

- Bootstrap confidence intervals are computed as in equation (2.6), using
  the inverse of Hhat from equation (2.8).

- In Hhat, v_P,n is implemented as se_mplp and v_n is implemented as se_us,
  the local-polynomial standard error for T_rd^*, including Q.plus and
  Q.minus, but without the mPLP bias-correction component.

- The prepivoting quantile levels use Hhat^{-1}(u) =
  Phi(mhat * Phi^{-1}(u)). With a Gaussian wild bootstrap, this agrees with
  the paper's result that the prepivoted interval is exactly the corresponding
  mPLP bias-corrected interval when the bootstrap statistic is conditionally
  Gaussian.

- The normal wild-bootstrap confidence interval is computed by actual
  fixed-regressor wild-bootstrap simulation, not by plugging in the analytic
  Gaussian quantiles. As boot.B grows, it should approach the corresponding
  mPLP bias-corrected interval up to Monte Carlo error.


New returned values when boot.ci = TRUE
---------------------------------------

- ci_mplp_boot
- boot.B
- boot.success
- boot.probs
- boot.mhat
- boot.wild
- se_us
- boot.t, when boot.keep = TRUE

