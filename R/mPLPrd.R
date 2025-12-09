######################################################################################
## mPLPrd: an R package for efficient RDD inference
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 1.01 (9.12.2025)
######################################################################################


mPLPrd <- function(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL, alpha = NULL, fast = NULL, g.loo = NULL) {
  
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
  
  
  if (!is.null(kersel) && is.character(kersel) && length(kersel) == 1 && !(kersel %in% c("uni", "uniform", "tri", "triangular", "epa", "epanechnikov"))){
    stop("Kernel function incorrectly specified.\n")
  } else if (is.null(kersel)){
    kersel <- "tri"
  }
  
  if (!is.null(res) && is.character(res) && length(res) == 1 && !(res %in% c("loo","cct-hc0","cct-hc1","cct-hc2","cct-hc3"))) {
    stop("Residuals incorrectly specified.\n")
  } else if (is.null(res)) {
    res <- "cct-hc3"
  }
  
  
  if (is.null(alpha))  {
    alpha  <- 0.05
  } else if (alpha <= 0 | alpha >=1) {
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
  
  
  ##########################  Initialization  ########################## 
  
  # Sample size
  n <- length(x)
  
  # Specify variables at each side of the cutoff
  xp <- x[x>=c]; yp <- y[x>=c]
  xm <- x[x<c];  ym <- y[x<c] 
  
  
  # Bw regularization   
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-c))[bwcheck]
    h        <- max(h, bw.min)
  }
  
  
  ##########################  ATE and Bias estimation    ########################## 
  
  # ATE estimation
  ghatp   <- LP.est(y=yp,x=xp,x0=c,h=h,p=p,v=0,kersel=kersel)
  ghatm   <- LP.est(y=ym,x=xm,x0=c,h=h,p=p,v=0,kersel=kersel)
  tau.hat <- ghatp - ghatm
  
  # Initialization of bias estimation
  if (res == "loo") {U <- (x-c)/(2*h)} else {U <- (x-c)/(h)}    
  U.Z              <- abs(U) < 1 
  nonzero_indices  <- which(U.Z != 0) 
  ghat_vec <- numeric(n)
  if (g.loo == TRUE) ghat_vec.loo <- numeric(n)
  if (fast == FALSE) C.xi <- numeric(n)
  
  # Loop for bias estimation
  for (ix in seq_along(nonzero_indices)) {
    i <- nonzero_indices[ix]
    x0i <- x[i]
    if (x0i >= c) {
      
      ghat_vec[i]     <- LP.est(y=yp,x=xp,x0=x0i,h=h,p=p,v=0,kersel=kersel)
      
      if (fast == FALSE) C.xi[i]  <- LP.est(y=((xp-x0i)/h)^(p+1),x=xp,x0=x0i,h=h,p=p,v=0,kersel=kersel)
      
      if (g.loo == TRUE) {
        
        xp.loo <- x*(x>=c);             yp.loo <- y*(x>=c)
        xp.loo <- xp.loo[-i];           yp.loo <- yp.loo[-i]
        xp.loo <- xp.loo[xp.loo != 0];  yp.loo <- yp.loo[yp.loo != 0]
        
        ghat_vec.loo[i] <- LP.est(y=yp.loo,x=xp.loo,x0=x0i,h=h,p=p,v=0,kersel=kersel)
        
      }
      
    } else if (x0i < c) {

      ghat_vec[i]     <- LP.est(y=ym,x=xm,x0=x0i,h=h,p=p,v=0,kersel=kersel)
      
      if (fast == FALSE) C.xi[i]  <- LP.est(y=((xm-x0i)/h)^(p+1),x=xm,x0=x0i,h=h,p=p,v=0,kersel=kersel)
      
      if (g.loo == TRUE) {
        
        xm.loo <- x*(x<c);              ym.loo <- y*(x<c)
        xm.loo <- xm.loo[-i];           ym.loo <- ym.loo[-i]
        xm.loo <- xm.loo[xm.loo != 0];  ym.loo <- ym.loo[ym.loo != 0]
        
        ghat_vec.loo[i] <- LP.est(y=ym.loo,x=xm.loo,x0=x0i,h=h,p=p,v=0,kersel=kersel)
        
      }
    }
  }
  
  # Estimation of Q at each side of the cutoff
  if (fast == FALSE) {
    Cp.LP <- LP.est(y=C.xi[x>=c],x=xp,x0=c,h=h,p=p,v=0,kersel=kersel)
    Cm.LP <- LP.est(y=C.xi[x<c],x=xm,x0=c,h=h,p=p,v=0,kersel=kersel)
    Cp    <- LP.est(y=((xp-c)/h)^(p+1),x=xp,x0=c,h=h,p=p,v=0,kersel=kersel)
    Cm    <- LP.est(y=((xm-c)/h)^(p+1),x=xm,x0=c,h=h,p=p,v=0,kersel=kersel)
    Qp    <- Cp/Cp.LP;  Qm    <- Cm/Cm.LP
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
    residualsp_gs <- ghat_vec.loo[x>=c] - c(ghatp)
    residualsm_gs <- ghat_vec.loo[x<c] - c(ghatm)
  } else if (g.loo == FALSE) {
    residualsp_gs <- ghat_vec[x>=c] - c(ghatp)
    residualsm_gs <- ghat_vec[x<c] - c(ghatm)
  }
  
  Bhatp_lp      <- LP.est(y=residualsp_gs,x=xp,x0=c,h=h,p=p,v=0,kersel=kersel)
  Bhatm_lp      <- LP.est(y=residualsm_gs,x=xm,x0=c,h=h,p=p,v=0,kersel=kersel)
  Bhat_mlp      <- Bhatp_lp*Qp - Bhatm_lp*Qm
  
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
    
    epshatp <- epshat[x>=c];   epshatm <- epshat[x< c]
  } else if (res == "cct-hc0" | res == "cct-hc1" | res == "cct-hc2" | res == "cct-hc3") {
    rp.pp1 <- matrix(NA, nrow = length(xp), ncol = (p + 2))
    rm.pp1 <- matrix(NA, nrow = length(xm), ncol = (p + 2))
    for (ip in 1:(p+2)) {
      rp.pp1[, ip] <- (xp-c)^(ip-1)
      rm.pp1[, ip] <- (xm-c)^(ip-1)
    }
    rp.p <- rp.pp1[,1:(p+1)]
    rm.p <- rm.pp1[,1:(p+1)]
    Kp.X    <- (K(u = (xp-c)/h, kersel = kersel)/h)
    Km.X    <- (K(u = (xm-c)/h, kersel = kersel)/h)
    invGp.q <- qrXXinv((sqrt(Kp.X)*rp.pp1))
    invGm.q <- qrXXinv((sqrt(Km.X)*rm.pp1))
    betap.q <- invGp.q%*%crossprod(rp.pp1*Kp.X,yp)
    betam.q <- invGm.q%*%crossprod(rm.pp1*Km.X,ym)
    
    if (res == "cct-hc0") {
      epshatp <- (yp - rp.pp1%*%betap.q)
      epshatm <- (ym - rm.pp1%*%betam.q)
      
      print(epshatp)
    } else {
      Qp.vec  <- rowSums((rp.pp1 %*% invGp.q) * (rp.pp1 * Kp.X))
      Qm.vec  <- rowSums((rm.pp1 %*% invGm.q) * (rm.pp1 * Km.X))
      if (res == "cct-hc1") {
        epshatp <- (yp - rp.pp1%*%betap.q)/(((length(yp) - p+1)/length(yp))^(0.5))
        epshatm <- (ym - rm.pp1%*%betam.q)/(((length(ym) - p+1)/length(ym))^(0.5))
      } else if (res == "cct-hc2") {
        epshatp <- (yp - rp.pp1%*%betap.q)/((1-Qp.vec)^(0.5))
        epshatm <- (ym - rm.pp1%*%betam.q)/((1-Qm.vec)^(0.5))
      } else if (res == "cct-hc3") {
        epshatp <- (yp - rp.pp1%*%betap.q)/((1-Qp.vec))
        epshatm <- (ym - rm.pp1%*%betam.q)/((1-Qm.vec))
      }
    }
  }
  
  # mPLP Standard errors
  if (fast == 0) {
    sep_mplp <- vmPLP(x=xp,x0=c,p=p,h=h,kersel=kersel,epshat = epshatp, Q=Qp)/sqrt(length(xp)*h) 
    sem_mplp <- vmPLP(x=xm,x0=c,p=p,h=h,kersel=kersel,epshat = epshatm, Q=Qm)/sqrt(length(xm)*h) 
  } else if (fast == 1) {
    if (kersel == "tri") {
      Kratio <- 0.84
    } else if (kersel == "epa") {
      Kratio <- 0.83 
    } else if (kersel == "uni") {
      Kratio <- 0.86 
    } 
    sep_mplp <-  Kratio*(CCT.se(x=xp,x0=c,h=h,p=p,kersel=kersel,epshat = epshatp)$se.rb)
    sem_mplp <-  Kratio*(CCT.se(x=xm,x0=c,h=h,p=p,kersel=kersel,epshat = epshatm)$se.rb)
  }
  
  se_mplp  <- (sqrt(sep_mplp^2 + sem_mplp^2))
  
  
  # Confidence Intervals
  ci_mplp <- c(tau.hat_mlpbc - qnorm(1-alpha/2)*se_mplp, tau.hat_mlpbc - qnorm(alpha/2)*se_mplp)

  
  # Output
  out <- list(n=n, h=h, Q.plus=Qp, Q.minus=Qm, tau.hat = tau.hat, tau.hat_mlpbc = tau.hat_mlpbc, 
              se_mplp=se_mplp, ci_mplp=ci_mplp, epshat.plus=epshatp, epshat.minus=epshatm)
  
  return(out)
  
}
