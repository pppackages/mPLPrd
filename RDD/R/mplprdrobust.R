mplprdrobust <- function(y, x, eval, h, p=NULL, kersel=NULL, res=NULL, Eg.loo = NULL, alpha = NULL, fast = NULL) {

   
  # The function mplprobust performs robust bias correction for a nonparametric
  #   regression via the modified prepivoted approach (mPLP*) of
  #   Cavaliere, Goncalves, Nielsen & Zanelli (2025)
  
  # y       <- dependent variable
  # x       <- regressor
  # eval    <- evaluation point
  # h       <- bandwidth
  # p       <- polynomial order      (default = 1)
  # kersel  <- kernel function used  (default = "tri")
  # res     <- residuals class used  (default = "cct")
  # alpha   <- significance level    (default = 0.05)
  
  if (is.null(p))         p      <- 1
  if (is.null(kersel))    kersel <- "tri"
  if (is.null(res))       res    <- "cct"
  if (is.null(Eg.loo))    Eg.loo <- FALSE
  if (is.null(alpha))     alpha  <- 0.05
  if (is.null(fast))      fast   <- FALSE
  
  deriv <- 0
  
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-eval))[bwcheck]
    h        <- max(h, bw.min)
  }
  
  # Sample size
  n <- length(x)
  
  # Specify regressors at the right and left of the cutoff
  xp <- x[x>=eval]
  yp <- y[x>=eval]
  xm <- x[x<eval]
  ym <- y[x<eval] 
  
  # Estimate tauh
  ghatp   <- LP.est(y=yp,x=xp,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
  ghatm   <- LP.est(y=ym,x=xm,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
  tau.hat <- ghatp - ghatm
  
  U <- (x-eval)/(h)
  if (res == "loo") {U <- (x-eval)/(2*h)}    
  if (Eg.loo == TRUE) ghat_vec.loo <- numeric(length(x))
  if (fast == FALSE) C.xi <- numeric(length(x))
  
  U.Z              <- abs(U) < 1 
  nonzero_indices  <- which(U.Z != 0) 
  ghat_vec         <- numeric(length(x))
  
  for (ix in seq_along(nonzero_indices)) {
    i <- nonzero_indices[ix]
    x0i <- x[i]
    if (x0i >= eval) {
      
      ghat_vec[i]     <- LP.est(y=yp,x=xp,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
      
      if (fast == FALSE)   C.xi[i]  <- LP.est(y=((xp-x0i)/h)^(p+1),x=xp,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
      
      if (Eg.loo == TRUE) {
        
        xp.loo <- x*(x>=eval)
        xp.loo <- xp.loo[-i]
        xp.loo <- xp.loo[xp.loo != 0]
        
        yp.loo <- y*(x>=eval)
        yp.loo <- yp.loo[-i]
        yp.loo <- yp.loo[yp.loo != 0]
        
        ghat_vec.loo[i] <- LP.est(y=yp.loo,x=xp.loo,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
        
      }
      
    } else if (x0i < eval) {

      ghat_vec[i]     <- LP.est(y=ym,x=xm,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
      
      if (fast == FALSE) C.xi[i]  <- LP.est(y=((xm-x0i)/h)^(p+1),x=xm,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
      
      if (Eg.loo == TRUE) {
        
        xm.loo <- x*(x<eval)
        xm.loo <- xm.loo[-i]
        xm.loo <- xm.loo[xm.loo != 0]
        
        ym.loo <- y*(x<eval)
        ym.loo <- ym.loo[-i]
        ym.loo <- ym.loo[ym.loo != 0]
        
        ghat_vec.loo[i] <- LP.est(y=ym.loo,x=xm.loo,x0=x0i,h=h,p=p,v=deriv,kersel=kersel)
        
      }
    }
  }
  
  # Correction
  if (fast == FALSE) {
    Cp.LP <- LP.est(y=C.xi[x>=eval],x=xp,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
    Cm.LP <- LP.est(y=C.xi[x<eval],x=xm,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
    Cp    <- LP.est(y=((xp-eval)/h)^(p+1),x=xp,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
    Cm    <- LP.est(y=((xm-eval)/h)^(p+1),x=xm,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
    Qp    <- Cp/Cp.LP
    Qm    <- Cm/Cm.LP
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
  
  # LP* bias correction
  if (Eg.loo == TRUE) {
    residualsp_gs <- ghat_vec.loo[x>=eval] - c(ghatp)
    residualsm_gs <- ghat_vec.loo[x<eval] - c(ghatm)
  } else if (Eg.loo == FALSE) {
    residualsp_gs <- ghat_vec[x>=eval] - c(ghatp)
    residualsm_gs <- ghat_vec[x<eval] - c(ghatm)
  }
  
  Bhatp_lp      <- LP.est(y=residualsp_gs,x=xp,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
  Bhatm_lp      <- LP.est(y=residualsm_gs,x=xm,x0=eval,h=h,p=p,v=deriv,kersel=kersel)
  Bhat_mlp      <- Bhatp_lp*Qp - Bhatm_lp*Qm
  tau.hat_mlpbc <- tau.hat - Bhat_mlp
  
  # Evaluate residuals
  if (res == "loo") {
    epshat <- y - ghat_vec.loo
    epshatp <- epshat[x>=eval]
    epshatm <- epshat[x< eval]
  } else if (res == "cct") {
    rp.pp1 <- matrix(NA, nrow = length(xp), ncol = (p + 2))
    rm.pp1 <- matrix(NA, nrow = length(xm), ncol = (p + 2))
    for (ip in 1:(p+2)) {
      rp.pp1[, ip] <- (xp-eval)^(ip-1)
      rm.pp1[, ip] <- (xm-eval)^(ip-1)
      
    }
    rp.p <- rp.pp1[,1:(p+1)]
    rm.p <- rm.pp1[,1:(p+1)]
    
    Kp.X <- (K(u = (xp-eval)/h, kersel = kersel)/h)
    Km.X <- (K(u = (xm-eval)/h, kersel = kersel)/h)
    invGp.q  <- qrXXinv((sqrt(Kp.X)*rp.pp1))
    invGm.q  <- qrXXinv((sqrt(Km.X)*rm.pp1))
    betap.q <- invGp.q%*%crossprod(rp.pp1*Kp.X,yp)
    betam.q <- invGm.q%*%crossprod(rm.pp1*Km.X,ym)
    Qp.vec <- rowSums((rp.pp1 %*% invGp.q) * (rp.pp1 * Kp.X))
    Qm.vec <- rowSums((rm.pp1 %*% invGm.q) * (rm.pp1 * Km.X))
    epshatp <- (yp - rp.pp1%*%betap.q)/((1-Qp.vec))
    epshatm <- (ym - rm.pp1%*%betam.q)/((1-Qm.vec))
  }
  
  # mPLP Standard error and confidence intervals
  if (fast == 0) {
    
    sep_mplp <- vmPLP(x=xp,x0=eval,p=p,h=h,kersel=kersel,epshat = epshatp, Q=Qp)/sqrt(length(xp)*h) 
    sem_mplp <- vmPLP(x=xm,x0=eval,p=p,h=h,kersel=kersel,epshat = epshatm, Q=Qm)/sqrt(length(xm)*h) 
    
  } else if (fast == 1) {
    
    if (kersel == "tri") {
      Kratio <- 0.84
    } else if (kersel == "epa") {
      Kratio <- 0.83 
    } else if (kersel == "uni") {
      Kratio <- 0.86 
    } 
    
    sep_mplp <-  Kratio*(CCT.se(x=xp,x0=eval,h=h,p=p,kersel=kersel,epshat = epshatp)$se.rb)
    sem_mplp <-  Kratio*(CCT.se(x=xm,x0=eval,h=h,p=p,kersel=kersel,epshat = epshatm)$se.rb)
    
  }
  
  se_mplp  <- (sqrt(sep_mplp^2 + sem_mplp^2))
  ci_mplp <- c(tau.hat_mlpbc - qnorm(1-alpha/2)*se_mplp, tau.hat_mlpbc - qnorm(alpha/2)*se_mplp)

  out <- list(n=n, h=h, Q.plus=Qp, Q.minus=Qm, tau.hat = tau.hat, tau.hat_mlpbc = tau.hat_mlpbc, 
              se_mplp=se_mplp, ci_mplp=ci_mplp, epshat.plus=epshatp, epshat.minus=epshatm)
  
  return(out)
  
}


