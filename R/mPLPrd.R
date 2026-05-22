######################################################################################
## mPLPrd: an R package for efficient RDD inference
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 1.04 (21.05.2026)
######################################################################################


mPLPrd <- function(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
                   alpha = NULL, fast = NULL, g.loo = NULL,
                   boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                   boot.wild = NULL, boot.keep = NULL) {

  # y       <- dependent variable
  # x       <- regressor
  # h       <- bandwidth
  # c       <- cutoff
  # p       <- polynomial order                 (default = 1)
  # kersel  <- kernel function used             (default = "tri")
  # res     <- residuals class used             (default = "cct-hc3")
  # alpha   <- significance level               (default = 0.05)
  # fast    <- fast std.error computation       (default = FALSE)
  # g.loo   <- leave-one-out Bias estimation    (default = FALSE)
  # boot.ci <- fixed-regressor wild bootstrap CI (default = FALSE)
  # boot.B  <- bootstrap repetitions             (default = 999 if boot.ci = TRUE)


  ############################# Error Checks #############################

  if (!is.numeric(h) || h <= 0) {
    stop("'h' must be a positive number")
  }


  if (!is.numeric(c) & !is.null(c)) {
    stop("'c' must be numeric")
  } else if (is.null(c)) {
    c <- 0
  }


  if (length(p) == 0) {
    p  <- 1
  } else if (length(p) > 1) {
    stop("Polynomial order p incorrectly specified.\n")
  }


  if (!is.null(kersel) && is.character(kersel) && length(kersel) == 1 && !(kersel %in% c("uni", "uniform", "tri", "triangular", "epa", "epanechnikov"))) {
    stop("Kernel function incorrectly specified.\n")
  } else if (is.null(kersel)) {
    kersel <- "tri"
  }
  kersel <- .canonical_kernel(kersel)

  if (!is.null(res) && is.character(res) && length(res) == 1 && !(res %in% c("loo", "cct-hc0", "cct-hc1", "cct-hc2", "cct-hc3"))) {
    stop("Residuals incorrectly specified.\n")
  } else if (is.null(res)) {
    res <- "cct-hc3"
  }


  if (is.null(alpha))  {
    alpha  <- 0.05
  } else if (alpha <= 0 | alpha >= 1) {
    stop("significance level incorrectly specified.\n")
  }


  if (is.null(fast)) {
    fast <- FALSE
  } else if (!is.logical(fast) || length(fast) > 1) {
    stop("'fast' must be a single TRUE/FALSE")
  }


  if (is.null(g.loo)) {
    g.loo <- FALSE
  } else if (!is.logical(g.loo) || length(g.loo) > 1) {
    stop("'g.loo' must be a single TRUE/FALSE")
  }


  if (is.null(boot.ci)) {
    boot.ci <- FALSE
  } else if (!is.logical(boot.ci) || length(boot.ci) > 1) {
    stop("'boot.ci' must be a single TRUE/FALSE")
  }


  if (is.null(boot.B)) {
    boot.B <- if (boot.ci) 999 else 0
  } else if (!is.numeric(boot.B) || length(boot.B) > 1 || boot.B < 1 || boot.B != floor(boot.B)) {
    stop("'boot.B' must be a positive integer")
  }
  boot.B <- as.integer(boot.B)


  if (!is.null(boot.seed) && (!is.numeric(boot.seed) || length(boot.seed) > 1)) {
    stop("'boot.seed' must be a single numeric value")
  }


  if (is.null(boot.wild)) {
    boot.wild <- "rademacher"
  } else if (!is.character(boot.wild) || length(boot.wild) > 1 ||
             !(boot.wild %in% c("normal", "rademacher"))) {
    stop("'boot.wild' must be 'normal' or 'rademacher'")
  }


  if (is.null(boot.keep)) {
    boot.keep <- TRUE
  } else if (!is.logical(boot.keep) || length(boot.keep) > 1) {
    stop("'boot.keep' must be a single TRUE/FALSE")
  }


  ##########################  Initialization  ##########################

  # Sample size
  n <- length(x)

  # Specify variables at each side of the cutoff
  indp <- x >= c
  indm <- x < c
  xp <- x[indp]
  yp <- y[indp]
  xm <- x[indm]
  ym <- y[indm]


  # Bw regularization
  h <- .regularize_h(x = x, x0 = c, h = h)


  ##########################  ATE and Bias estimation    ##########################

  # ATE estimation
  ghatp   <- LP.est(y = yp, x = xp, x0 = c, h = h, p = p, v = 0, kersel = kersel)
  ghatm   <- LP.est(y = ym, x = xm, x0 = c, h = h, p = p, v = 0, kersel = kersel)
  tau.hat <- ghatp - ghatm

  # Initialization of bias estimation
  # The mPLP bias correction uses auxiliary local fits over the 2h window.
  U <- (x - c) / (2 * h)
  U.Z <- abs(U) < 1
  nonzero_indices <- which(U.Z != 0)
  ghat_vec <- numeric(n)
  if (g.loo == TRUE) ghat_vec.loo <- numeric(n)
  if (fast == FALSE) C.xi <- numeric(n)

  # Loop for bias estimation
  for (ix in seq_along(nonzero_indices)) {
    i <- nonzero_indices[ix]
    x0i <- x[i]
    if (x0i >= c) {

      ghat_vec[i] <- LP.est(y = yp, x = xp, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      if (fast == FALSE) C.xi[i] <- LP.est(y = ((xp - x0i) / h)^(p + 1), x = xp, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      if (g.loo == TRUE) {

        xp.loo <- x * (x >= c)
        yp.loo <- y * (x >= c)
        xp.loo <- xp.loo[-i]
        yp.loo <- yp.loo[-i]
        xp.loo <- xp.loo[xp.loo != 0]
        yp.loo <- yp.loo[yp.loo != 0]

        ghat_vec.loo[i] <- LP.est(y = yp.loo, x = xp.loo, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      }

    } else if (x0i < c) {

      ghat_vec[i] <- LP.est(y = ym, x = xm, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      if (fast == FALSE) C.xi[i] <- LP.est(y = ((xm - x0i) / h)^(p + 1), x = xm, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      if (g.loo == TRUE) {

        xm.loo <- x * (x < c)
        ym.loo <- y * (x < c)
        xm.loo <- xm.loo[-i]
        ym.loo <- ym.loo[-i]
        xm.loo <- xm.loo[xm.loo != 0]
        ym.loo <- ym.loo[ym.loo != 0]

        ghat_vec.loo[i] <- LP.est(y = ym.loo, x = xm.loo, x0 = x0i, h = h, p = p, v = 0, kersel = kersel)

      }
    }
  }

  # Estimation of Q at each side of the cutoff
  if (fast == FALSE) {
    Cp.LP <- LP.est(y = C.xi[indp], x = xp, x0 = c, h = h, p = p, v = 0, kersel = kersel)
    Cm.LP <- LP.est(y = C.xi[indm], x = xm, x0 = c, h = h, p = p, v = 0, kersel = kersel)
    Cp    <- LP.est(y = ((xp - c) / h)^(p + 1), x = xp, x0 = c, h = h, p = p, v = 0, kersel = kersel)
    Cm    <- LP.est(y = ((xm - c) / h)^(p + 1), x = xm, x0 = c, h = h, p = p, v = 0, kersel = kersel)
    Qp    <- Cp / Cp.LP
    Qm    <- Cm / Cm.LP
  } else if (fast == TRUE) {
    if (kersel == "tri") {
      Qp <- 1.4082
    } else if (kersel == "epa") {
      Qp <- 1.3571
    } else if (kersel == "uni") {
      Qp <- 1.2
    }
    Qm <- Qp
  }

  # Bias Estimation
  if (g.loo == TRUE) {
    residualsp_gs <- ghat_vec.loo[indp] - c(ghatp)
    residualsm_gs <- ghat_vec.loo[indm] - c(ghatm)
  } else if (g.loo == FALSE) {
    residualsp_gs <- ghat_vec[indp] - c(ghatp)
    residualsm_gs <- ghat_vec[indm] - c(ghatm)
  }

  Bhatp_lp      <- LP.est(y = residualsp_gs, x = xp, x0 = c, h = h, p = p, v = 0, kersel = kersel)
  Bhatm_lp      <- LP.est(y = residualsm_gs, x = xm, x0 = c, h = h, p = p, v = 0, kersel = kersel)
  Bhat_mlp      <- Bhatp_lp * Qp - Bhatm_lp * Qm

  # De-biased ATE
  tau.hat_mlpbc <- tau.hat - Bhat_mlp


  #####################  Standard errors and Confidence Intervals  #####################

  # Residuals
  if (res == "loo") {
    if (g.loo == TRUE) {
      epshat <- y - ghat_vec.loo
    } else {
      epshat <- y - ghat_vec
    }

    epshatp <- epshat[indp]
    epshatm <- epshat[indm]
  } else if (res == "cct-hc0" | res == "cct-hc1" | res == "cct-hc2" | res == "cct-hc3") {
    rp.pp1 <- .poly_basis(xp - c, p + 1)
    rm.pp1 <- .poly_basis(xm - c, p + 1)
    Kp.X <- K(u = (xp - c) / h, kersel = kersel) / h
    Km.X <- K(u = (xm - c) / h, kersel = kersel) / h
    invGp.q <- qrXXinv(sqrt(Kp.X) * rp.pp1)
    invGm.q <- qrXXinv(sqrt(Km.X) * rm.pp1)
    betap.q <- invGp.q %*% crossprod(rp.pp1 * Kp.X, yp)
    betam.q <- invGm.q %*% crossprod(rm.pp1 * Km.X, ym)

    if (res == "cct-hc0") {
      epshatp <- yp - rp.pp1 %*% betap.q
      epshatm <- ym - rm.pp1 %*% betam.q
    } else {
      Qp.vec <- rowSums((rp.pp1 %*% invGp.q) * (rp.pp1 * Kp.X))
      Qm.vec <- rowSums((rm.pp1 %*% invGm.q) * (rm.pp1 * Km.X))
      if (res == "cct-hc1") {
        epshatp <- (yp - rp.pp1 %*% betap.q) / (((length(yp) - p + 1) / length(yp))^(0.5))
        epshatm <- (ym - rm.pp1 %*% betam.q) / (((length(ym) - p + 1) / length(ym))^(0.5))
      } else if (res == "cct-hc2") {
        epshatp <- (yp - rp.pp1 %*% betap.q) / ((1 - Qp.vec)^(0.5))
        epshatm <- (ym - rm.pp1 %*% betam.q) / ((1 - Qm.vec)^(0.5))
      } else if (res == "cct-hc3") {
        epshatp <- (yp - rp.pp1 %*% betap.q) / (1 - Qp.vec)
        epshatm <- (ym - rm.pp1 %*% betam.q) / (1 - Qm.vec)
      }
    }
  }

  # mPLP Standard errors
  if (fast == 0) {
    sep_mplp <- vmPLP(x = xp, x0 = c, p = p, h = h, kersel = kersel, epshat = epshatp, Q = Qp) / sqrt(length(xp) * h)
    sem_mplp <- vmPLP(x = xm, x0 = c, p = p, h = h, kersel = kersel, epshat = epshatm, Q = Qm) / sqrt(length(xm) * h)
  } else if (fast == 1) {
    if (kersel == "tri") {
      Kratio <- 0.84
    } else if (kersel == "epa") {
      Kratio <- 0.83
    } else if (kersel == "uni") {
      Kratio <- 0.86
    }
    sep_mplp <- Kratio * (CCT.se(x = xp, x0 = c, h = h, p = p, kersel = kersel, epshat = epshatp)$se.rb)
    sem_mplp <- Kratio * (CCT.se(x = xm, x0 = c, h = h, p = p, kersel = kersel, epshat = epshatm)$se.rb)
  }

  se_mplp <- sqrt(sep_mplp^2 + sem_mplp^2)


  # Confidence Intervals
  ci_mplp <- c(tau.hat_mlpbc - qnorm(1 - alpha / 2) * se_mplp, tau.hat_mlpbc - qnorm(alpha / 2) * se_mplp)


  #####################  Fixed-regressor wild bootstrap CI  #####################

  if (boot.ci == TRUE) {
    if (!is.null(boot.seed)) set.seed(as.integer(boot.seed))

    epshat <- numeric(n)
    epshat[indp] <- as.numeric(epshatp)
    epshat[indm] <- as.numeric(epshatm)

    sep_us <- Qp * LP.se(x = xp, x0 = c, h = h, p = p, kersel = kersel, epshat = epshatp)
    sem_us <- Qm * LP.se(x = xm, x0 = c, h = h, p = p, kersel = kersel, epshat = epshatm)
    se_us <- sqrt(sep_us^2 + sem_us^2)

    mhat <- as.numeric(se_mplp / se_us)
    if (!is.finite(mhat) || mhat <= 0) {
      stop("Bootstrap CI cannot be computed because the non-bias-corrected standard error is not positive.")
    }

    H.inv <- function(u) pnorm(mhat * qnorm(u))
    boot.probs <- H.inv(c(alpha / 2, 1 - alpha / 2))

    h.boot.p <- .regularize_h(x = xp, x0 = c, h = h)
    h.boot.m <- .regularize_h(x = xm, x0 = c, h = h)
    boot.local <- (indp & abs(x - c) <= h.boot.p) |
      (indm & abs(x - c) <= h.boot.m)
    if (!any(boot.local & indp) || !any(boot.local & indm)) {
      stop("Bootstrap CI requires observations on both sides of the cutoff inside the bandwidth.")
    }

    ghat.boot.mean <- numeric(n)
    boot.local.indices <- which(boot.local)
    for (ix in seq_along(boot.local.indices)) {
      i <- boot.local.indices[ix]
      if (x[i] >= c) {
        ghat.boot.mean[i] <- LP.est(y = yp, x = xp, x0 = x[i], h = h, p = p, v = 0, kersel = kersel)
      } else {
        ghat.boot.mean[i] <- LP.est(y = ym, x = xm, x0 = x[i], h = h, p = p, v = 0, kersel = kersel)
      }
    }

    boot.t <- rep(NA_real_, boot.B)

    for (bb in seq_len(boot.B)) {
      e.boot <- if (boot.wild == "normal") {
        rnorm(length(boot.local.indices))
      } else {
        sample(c(-1, 1), length(boot.local.indices), replace = TRUE)
      }

      y.boot <- y
      y.boot[boot.local.indices] <- ghat.boot.mean[boot.local.indices] +
        epshat[boot.local.indices] * e.boot

      ghatp.boot <- LP.est(y = y.boot[indp], x = xp, x0 = c, h = h, p = p, v = 0, kersel = kersel)
      ghatm.boot <- LP.est(y = y.boot[indm], x = xm, x0 = c, h = h, p = p, v = 0, kersel = kersel)

      boot.t[bb] <- as.numeric(Qp * (ghatp.boot - ghatp) - Qm * (ghatm.boot - ghatm))
    }

    boot.t.ok <- boot.t[is.finite(boot.t)]
    boot.success <- length(boot.t.ok)
    if (boot.success < 2) {
      warning("Fewer than two successful bootstrap repetitions; bootstrap CI is NA.")
      ci_mplp_boot <- c(NA_real_, NA_real_)
    } else {
      q.boot <- as.numeric(quantile(boot.t.ok, probs = boot.probs, names = FALSE, type = 7))
      ci_mplp_boot <- c(as.numeric(tau.hat) - q.boot[2], as.numeric(tau.hat) - q.boot[1])
    }
  }

  
  # Output
  out <- list(n = n, h = h, Q.plus = Qp, Q.minus = Qm, tau.hat = tau.hat, tau.hat_mlpbc = tau.hat_mlpbc,
              se_mplp = se_mplp, ci_mplp = ci_mplp, epshat.plus = epshatp, epshat.minus = epshatm)

  if (boot.ci == TRUE) {
    out <- c(out, list(ci_mplp_boot = ci_mplp_boot,
                       boot.B = boot.B,
                       boot.success = boot.success,
                       boot.probs = boot.probs,
                       boot.mhat = mhat,
                       boot.wild = boot.wild,
                       se_us = se_us))
    if (boot.keep == TRUE) out$boot.t <- boot.t
  }
  
  return(out)

}
