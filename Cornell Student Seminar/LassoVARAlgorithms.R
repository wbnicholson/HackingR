# This script contains three versions of the same Lasso VAR algorithm with various optimizations.

# Requires the Rcpp and RcppArmadillo

library(Rcpp)
library(RcppArmadillo)

# First Lasso-VAR algorithm.  Consists of three for loops.
# Inner "kernel" of Lasso-VAR function

lasso <- function(Z, Y, gam, maxiter=1e3, tol, B=NULL,znorm2) {
  pk <- nrow(Z)
  T <- ncol(Z)
  k <- nrow(Y)
  if (is.null(B)) B <- matrix(0, nrow=k, ncol=pk)
  for (i in seq(maxiter)) {
    BOLD <- B
    for (l in seq(nrow(B))) {
      for (m in seq(ncol(B))) {
        r <- Y[l,] - B[l,-m]%*%Z[-m,] 
        B[l, m] <- ST(sum(r*Z[m,]), gam) / znorm2[m]
      }
    }
    if (max(abs(B-BOLD)/(1+BOLD)) < tol) { 
      return(B)
    }
  }
}

lasso2 <- function(Z, Y, gam, maxiter=1e3, tol, B=NULL,znorm2) {
  pk <- nrow(Z)
  T <- ncol(Z)
  k <- nrow(Y)
  if (is.null(B)) B <- matrix(0, nrow=k, ncol=pk)
  for (i in seq(maxiter)) {
    BOLD <- B
    for (l in seq(nrow(B))) {
      for (m in seq(ncol(B))) {
        r <- Y[l,] - B[l,-m]%*%Z[-m,] 
        B[l, m] <- ST1(sum(r*Z[m,]), gam) / znorm2[m]
      }
    }
    if (max(abs(B-BOLD)/(1+BOLD)) < tol) { 
      return(B)
    }
  }
}


# Algorithm used for a grid of lambda values with scaling
lassoVAR <- function (B, Z, Y, gamm, eps) 
{
    beta = list()
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    znorm2 <- rowSums(Z^2)
    for (i in 1:length(gamm)) {
        B[[i]] <- matrix(B[[i]][, 2:ncol(B[[i]])],nrow=k,ncol=k*p)
    }
    BOLD=B[[1]]
    for (i in 1:length(gamm)) {
        gam <- gamm[i]
        B <- lasso(Z, Y, gam,1e3,eps,B=BOLD,znorm2)
        nu <- c(apply(YOLD, 1, mean)) - B %*% apply(ZOLD, 1, 
            mean)
        beta[[i]] <- cbind(nu, B)
    }
    return(list(beta = beta, gamma = gamm))
}

lassoVAR2 <- function (B, Z, Y, gamm, eps) 
{
  beta = list()
  Y <- t(Y)
  YOLD <- Y
  ZOLD <- Z
  Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
  Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  znorm2 <- rowSums(Z^2)
  for (i in 1:length(gamm)) {
    B[[i]] <- matrix(B[[i]][, 2:ncol(B[[i]])],nrow=k,ncol=k*p)
  }
  BOLD=B[[1]]
  for (i in 1:length(gamm)) {
    gam <- gamm[i]
    B <- lasso2(Z, Y, gam,1e3,eps,B=BOLD,znorm2)
    nu <- c(apply(YOLD, 1, mean)) - B %*% apply(ZOLD, 1, 
                                                mean)
    beta[[i]] <- cbind(nu, B)
  }
  return(list(beta = beta, gamma = gamm))
}

# Lasso "kernel" which uses vectorized operations:


lassocorevec <- function(Z, Y, gam, maxiter=1e3, tol, B=NULL,znorm2) {
  pk <- nrow(Z)
  T <- ncol(Z)
  k <- nrow(Y)
  if (is.null(B)) B <- matrix(0, nrow=k, ncol=pk)
  for (i in seq(maxiter)) {
    BOLD <- B
      for (m in seq(ncol(B))) {
        r <- Y - B[,-m]%*%Z[-m,] # partial residual
        B[, m] <- ST3a(rowSums(r%*%Z[m,]), gam) / znorm2[m]
      }
    if (max(abs(B-BOLD)/(1+BOLD)) < tol) { # check if converged
      return(B)
    }
  }
}

# Full Lasso VAR with vectorized operations
lassoVARVec <- function (B, Z, Y, gamm, eps) 
{
    beta = list()
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    znorm2 <- rowSums(Z^2)
    for (i in 1:length(gamm)) {
        B[[i]] <- matrix(B[[i]][, 2:ncol(B[[i]])],nrow=k,ncol=k*p)
    }
    BOLD <- B[[1]]
    for (i in 1:length(gamm)) {
        gam <- gamm[i]
        B <- lassocorevec(Z, Y, gam,1e3,eps,B=BOLD,znorm2)
        nu <- c(apply(YOLD, 1, mean)) - B %*% apply(ZOLD, 1, 
            mean)
        beta[[i]] <- cbind(nu, B)
    }
    return(list(beta = beta, gamma = gamm))
}
sourceCpp("lassocoreserial.cpp")
# Lasso VAR function with C++ based kernel
lassoVARCPP <- function (B, Z, Y, gamm, eps) 
{
    beta = list()
    Y <- t(Y)

    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    znorm2 <- rowSums(Z^2)
    for (i in 1:length(gamm)) {
        B[[i]] <- matrix(B[[i]][, 2:ncol(B[[i]])],nrow=k,ncol=k*p)
    }
    beta <- gamloop(B,Y,Z,gamm,eps,as.matrix(YMean),as.matrix(ZMean),znorm2,B[[1]])

    return(list(beta = beta, gamma = gamm))
}

## Elementwise soft-thresholding
ST <- function(x,gam)
    {
        return(sign(x)*max(abs(x)-gam,0))

        }
# vector-wise soft-thresholding
ST3a <- function(x,gam)
    {

        return(sign(x)*pmax(abs(x)-rep(gam,length(x)),0))
        }


lassoVARVec <- function (B, Z, Y, gamm, eps) 
{
    beta = list()
    Y <- t(Y)
    YOLD <- Y
    ZOLD <- Z
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    znorm2 <- rowSums(Z^2)
    for (i in 1:length(gamm)) {
        B[[i]] <- matrix(B[[i]][, 2:ncol(B[[i]])],nrow=k,ncol=k*p)
    }
    BOLD <- B[[1]]
    for (i in 1:length(gamm)) {
        gam <- gamm[i]
        B <- lassocorevec(Z, Y, gam,1e3,eps,B=BOLD,znorm2)
        nu <- c(apply(YOLD, 1, mean)) - B %*% apply(ZOLD, 1, 
            mean)
        beta[[i]] <- cbind(nu, B)
    }
    return(list(beta = beta, gamma = gamm))
}
sourceCpp("lassocoreserial.cpp")
sourceCpp("lassocoreparallel.cpp")

# Lasso VAR function with parallel C++ based kernel
lassoVARCPPPar <- function (B, Z, Y, gamm, eps) 
{
    Y <- t(Y)
    YMean <- c(apply(Y, 1, mean))
    ZMean <- c(apply(Z, 1, mean))
    Y <- Y - c(apply(Y, 1, mean)) %*% t(c(rep(1, ncol(Y))))
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
    znorm2 <- rowSums(Z^2)
    BFOO <- B[, 2:ncol(B[, , 1]), 1]
    beta <- gamloopP(B[,2:ncol(B[,,1]),],Y,Z,gamm,eps,as.matrix(YMean),as.matrix(ZMean),znorm2,BFOO)
    return(list(beta = beta, gamma = gamm))
}

