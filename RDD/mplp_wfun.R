library(sandwich)

## -- Kernel function -- ## 

K <- function(u, kersel) {
  kersel <- kersel[1]  # force scalar
  
  result <- if (kersel == "epa") {
    ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
  } else if (kersel == "tri") {
    ifelse(abs(u) <= 1, (1 - abs(u)), 0)
  } else if (kersel == "uni") {
    ifelse(abs(u) <= 1, 0.5, 0)
  } else {
    warning("Unknown kernel type. Returning 0.")
    0
  }
  
  return(as.numeric(result))  
}

## -- polynomial basis function -- ## 

r_p <- function(p) {
  function(u) {
    # Ensure u is treated as a vector
    u <- as.vector(u)
    
    # Create the matrix: each row corresponds to a u value, each column u^j
    outer(u, 0:p, `^`)
  }
}

## -- Gamma matrix -- ##

Gamma.p.inv <- function(u,h,p,kersel) {
  
  X <- matrix(NA, nrow = length(u), ncol = p+1)
  for (i in 1:(p+1)) {
    X[,i] <- (u)^(i-1)
  }
  K.X <- K(u = u, kersel = kersel)/(h*length(u))
  nonzero_indices <- which(K.X != 0) 
  X <- X[nonzero_indices, ]
  K.X <- K.X[nonzero_indices]
  Gamma.inv <- qrXXinv((sqrt(K.X)*X)) 
}

## -- Local polynomial estimator -- ## 

LP.est <- function(y,x,x0,h,p,v,kersel) {
  
  n <- length(x)
  X <- matrix(NA, nrow = n, ncol = p+1)
  
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-x0))[bwcheck]
    h     <- max(h, bw.min)
  }
  
  for (i in 1:(p+1)) {
    X[,i] <- ((x-x0)/h)^(i-1)
  }
  K.X <- K(u = (x-x0)/h, kersel = kersel)/h
  nonzero_indices <- which(K.X != 0) 
  X <- X[nonzero_indices, ]
  y <- y[nonzero_indices]
  K.X <- K.X[nonzero_indices]
  Gamma.inv <- qrXXinv((sqrt(K.X)*X)) 
  
  
  e1 <- numeric(p+1)
  e1[v+1] <- 1
  ghat <- factorial(v)*e1%*%Gamma.inv%*%crossprod(X*K.X,y)/(h^v)
  return(ghat)
  
}


## -- Equivalent kernel of a local polynomial estimator evaluated at a single x_i -- ## 

wix0 <- function(x,xi,x0,h,p,v,kersel) {
  
  bwcheck <- 21
  if (!is.null(bwcheck)) {
    bw.min   <- sort(abs(x-x0))[bwcheck]
    h     <- max(h, bw.min)
  }
  
  e1 <- numeric(p+1)
  e1[v+1] <- 1
  n <- length(x)
  u <- (x-x0)/h
  X <- matrix(NA, nrow = n, ncol = p+1)
  for (i in 1:(p+1)) {
    X[,i] <- ((x-x0)/h)^(i-1)
  }
  K.X <- K(u = (x-x0)/h, kersel = kersel)/h
  nonzero_indices <- which(K.X != 0) 
  X <- X[nonzero_indices, ]
  K.X <- K.X[nonzero_indices]
  Gamma.inv <- qrXXinv((sqrt(K.X)*X))*n
  
  W <- numeric(length(xi))
  u_all <- (xi - x0) / h
  K_all <- K(u = u_all, kersel = kersel)
  nonzero_indices.2 <- which(K_all != 0)
  for (ii in nonzero_indices.2) {
    ui <- u_all[ii]
    kval <- K_all[ii]
    
    ru <- kval * ui^(0:p)
    W[ii] <- factorial(v) * e1 %*% Gamma.inv %*% ru
  }
  
  
  return(W)
}

## -- Equivalent kernel of a (mPLP based) debiased LP estimator evaluated at a single x_i -- ## 

wmPLP <- function(x,xi,x0,h,p,kersel, Q) {
  n <- length(x)
  A <- c(1+Q)*wix0(x = x,xi = xi, x0 = x0, h = h, p=p, v=0, kersel = kersel) 
  B <- numeric(n)
  
  U.Z = abs((x-x0)/h) <= 1
  nonzero_indices <- which(U.Z != 0) 
  for (jj in nonzero_indices) {
    
    BB <- c(Q) * wix0(x=x, xi=x[jj], x0=x0, h=h, p=p, v=0, kersel=kersel)*wix0(x=x, xi=xi, x0=x[jj], h=h, p=p, v=0, kersel=kersel)/(n*h)
    B[jj] <- BB
  } 
  B <- sum(B)
  
  A - B
  
}

## -- Variance estimator of a (mPLP-based) bias corrected test statistic -- ## 

vmPLP <- function(x,x0,h,p,kersel,epshat,Q) {
  
  n <- length(x)
  C <- numeric(n)
  U.Z <- (abs((x-x0)/(2*h)) <= 1) 
  nonzero_indices <- which(U.Z != 0) 
  for (ix in seq_along(nonzero_indices)) {
    
    ii <- nonzero_indices[ix]
    C[ii] <- ((wmPLP(x = x, xi = x[ii], x0 = x0, h=h, p=p, kersel=kersel, Q=Q)*epshat[ii])^2)/(n*h)
    
  } 
  C <- sum(C)
  
  sqrt(C)
}

## -- CCT standard errors -- ## 

CCT.se <- function(x,x0,h,p,kersel,epshat) {
  
  n <- length(x)
  r.pp1 <- matrix(NA, nrow = n, ncol = (p + 2))
  for (ip in 1:(p+2)) {
    r.pp1[, ip] <- (x-x0)^(ip-1)
  }
  r.p <- r.pp1[,1:(p+1)]
  
  e.p1    <- matrix(0,(p+2),1); e.p1[p+2]=1
  K.X <- K(u = (x-x0)/h, kersel = kersel)/h
  invG.q  <- qrXXinv((sqrt(K.X)*r.pp1))
  invG.p  <- qrXXinv((sqrt(K.X)*r.p))
  L <- crossprod(r.p*K.X,((x-x0)/h)^(p+1))
  
  Q.q <- t(t(r.p*K.X) - (h^(p+1))*(L%*%t(e.p1))%*%t(t(invG.q%*%t(r.pp1))*K.X))
  Q.p <- r.p*K.X
  
  MM.p <- crossprod(c(epshat)*Q.p)
  MM.q <- crossprod(c(epshat)*Q.q)
  MMM.p <- invG.p%*%MM.p%*%invG.p
  MMM.q <- invG.p%*%MM.q%*%invG.p
  
  
  M.p <- sqrt(MMM.p[1,1])
  M.q <- sqrt(MMM.q[1,1])
  
  ses <- list(se.us = M.p, se.rb=M.q)
  return(ses)
}

qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}


## -- CCT (HC3) residuals -- ## 


CCT.res = function(y,x,eval,h,p,kersel) {
  
  r.pp1 <- matrix(NA, nrow = n, ncol = (p + 2))
  for (ip in 1:(p+2)) {
    r.pp1[, ip] <- (x-eval)^(ip-1)
    
  }
  r.p <- r.pp1[,1:(p+1)]
  
  K.X <- K(u = (x-eval)/h, kersel = kersel)/h
  invG.q  <- qrXXinv((sqrt(K.X)*r.pp1))
  invG.p  <- qrXXinv((sqrt(K.X)*r.p))
  beta.q <- invG.q%*%crossprod(r.pp1*K.X,y)
  
  Q.res <- r.p%*%invG.p%*%t(r.p*K.X)
  Q.vec <- diag(Q.res)
  epshat <- (y - r.pp1%*%beta.q)/((1-Q.vec))
  
}



w_conv.int <- function(u, e1.Gamma.p.inv, e1.Gamma.p.inv.pop, p, kersel) {
  
  # function to integrate over r
  w_conv.inner <- function(r) {
    term1 <- as.numeric(e1.Gamma.p.inv %*% t(r_p(p)(r))) * K(r, kersel)
    term2 <- as.numeric(e1.Gamma.p.inv.pop %*% t(r_p(p)(u-r))) * K(u - r, kersel)
    term1 * term2
  }
  
  integrate(w_conv.inner, lower = -1, upper = 1)$value
}

make_w_conv <- function(e1.Gamma.p.inv, e1.Gamma.p.inv.pop, p, kersel) {
  function(u) {
    w_conv.int(u, e1.Gamma.p.inv, e1.Gamma.p.inv.pop, p, kersel)
  }
}



# w_conv.bnd <- function(u, e1.Gamma.p.inv, p, kersel) {
# 
#   e1    <- numeric(p+1)
#   e1[1] <- 1
# 
#   # function to integrate over r
#   w_conv.inner.bnd <- function(r) {
#     term1 <- as.numeric(e1.Gamma.p.inv %*% t(r_p(p)(r))) * K(r, kersel)
#     term2 <- as.numeric(e1%*%Gamma.p.pop.inv.fun(kersel=kersel, p=p)(r)%*%t(r_p(p)(u-r))) * K(u - r, kersel)
#     term1 * term2
#   }
# 
#   integrate(w_conv.inner.bnd, lower = 0, upper = 1)$value
# }

w_conv.bnd <- function(u, e1.Gamma.p.inv, p, kersel) {
  e1 <- numeric(p + 1); e1[1] <- 1
  
  # get the inverse-matrix function once
  Ginv_fun <- Gamma.p.pop.inv.fun(kersel = kersel, p = p)
  
  # integrand that handles vector r by computing elementwise
  w_conv.inner.bnd <- function(r_vec) {
    # r_vec may be scalar or a vector. We return numeric of same length.
    # compute term1 and term2 elementwise
    r_vec <- as.numeric(r_vec)
    # prepare polynomial evaluations:
    Rmat_r   <- r_p(p)(r_vec)        # length(r_vec) x (p+1)
    Rmat_ur  <- r_p(p)(u - r_vec)    # same dims
    
    # term1: e1.Gamma.p.inv %*% t(r_p(p)(r))  -> for each r: scalar
    # We compute row-wise multiplication: (Rmat_r %*% as.numeric(e1.Gamma.p.inv))
    # but e1.Gamma.p.inv is (p+1)-vector? In your code you use it as row vector: adapt accordingly.
    # Here I assume e1.Gamma.p.inv is a (1 x (p+1)) numeric vector (e1 %*% Gamma.inv)
    vec_term1 <- as.numeric(Rmat_r %*% as.numeric(e1.Gamma.p.inv)) * K(r_vec, kersel)
    
    # term2: need e1 %*% Gamma.p.pop.inv(scalar r) %*% t(r_p(p)(u-r))
    # For each r_i call Ginv_fun(r_i) (which returns a (p+1)x(p+1) matrix), then compute scalar
    vec_term2 <- vapply(seq_along(r_vec), function(ii) {
      r_i <- r_vec[ii]
      # get matrix Gamma^{-1} evaluated at s = r_i
      Ginv_mat <- Ginv_fun(r_i)  # <-- Ginv_fun expects scalar
      # r_p(u - r_i) is a row vector (1 x (p+1))
      rp_ur <- r_p(p)(u - r_i)   # will be 1 x (p+1)
      # compute e1 %*% Ginv_mat %*% t(rp_ur)  (e1 is 1 x (p+1))
      as.numeric(e1 %*% Ginv_mat %*% t(rp_ur)) * K(u - r_i, kersel)
    }, numeric(1))
    
    # final integrand values
    vec_term1 * vec_term2
  }
  
  # integrate from 0 to 1
  integrate(w_conv.inner.bnd, lower = 0, upper = 1)$value
}


make_w_conv_bnd <- function(e1.Gamma.p.inv, p, kersel) {
  function(u) {
    w_conv.bnd(u, e1.Gamma.p.inv, p, kersel)
  }
}


Gamma.p.pop <- function(kersel, p) {
  
  # Compute moments m_0 through m_{p+1}
  moments <- numeric(p + 2)
  
  for (k in 0:(2*p)) {
    moment_fun <- function(x) x^k * K(x,kersel)
    moments[k + 1] <- integrate(moment_fun, -1, 1)$value
  }
  
  # Build (p+1) x (p+1) matrix
  M <- matrix(0, nrow = p + 1, ncol = p + 1)
  
  for (i in 1:(p + 1)) {
    for (j in 1:(p + 1)) {
      M[i, j] <- moments[i + j - 1]   # i+j-1 corresponds to moment index (i+j-2)
    }
  }
  
  return(M)
}



Gamma.p.pop.inv.fun <- function(kersel, p) {
  
  # Returned function: evaluates inverse matrix at given s
  Minv_fun <- function(s) {
    
    # 1. Compute moments m_0 ... m_{p+1} over [-s, 1]
    moments <- numeric(p + 2)
    for (k in 0:(2*p)) {
      moment_fun     <- function(x) x^k * K(x,kersel)
      moments[k + 1] <- integrate(moment_fun, lower = -s , upper = 1)$value #integrate(moment_fun, lower = max(-s, -1), upper = 1)$value
    }
    
    # 2. Build (p+1)x(p+1) moment matrix
    M <- matrix(0, nrow = p + 1, ncol = p + 1)
    for (i in 1:(p + 1)) {
      for (j in 1:(p + 1)) {
        M[i, j] <- moments[i + j - 1]
      }
    }
    
    # 3. Invert and return
    return(solve(M))
  }
  
  return(Minv_fun)
}


fast.comb.rdrobust <- function(tauhat.mplp, tauhat.rbc, se.mplp, se.rbc,  kersel=NULL, alpha = NULL) {
  
  if (is.null(kersel))    kersel <- "tri"
  
  if (kersel == "tri") {
    lambda <- 1.4468
    covariance <- (0.7688*((se.rbc^2)))
  } else if (kersel == "uni") {
    lambda <- 0.7102
    covariance <- 0.5533*(se.rbc/2)
  } else if (kersel == "epa") {
    lambda <- 1.0811
    covariance <- 0.7127*(se.rbc/2)
  }
  
  
  tauhat.bc.comb <- lambda*tauhat.mplp + (1-lambda)*tauhat.rbc
  se.bc.comb <- sqrt( (lambda^2)*(se.mplp^2) + ((1-lambda)^2)*(se.rbc^2)  + 2*(1-lambda)*lambda*covariance )
  ci.comb <- c(tauhat.bc.comb - qnorm(1-alpha/2)*se.bc.comb, tauhat.bc.comb - qnorm(alpha/2)*se.bc.comb)
  
  out <- list(tauhat.bc.comb = tauhat.bc.comb, se.bc.comb=se.bc.comb, ci.comb=ci.comb)
  
  
}











