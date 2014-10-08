# converts a VAR of order p to a VAR of order 1 for simulation purposes
VarptoVar1 <- function (Ai, p, k) 
{
    Y <- matrix(0, nrow = k * p, ncol = 1)
    Y[1, ] <- 0
    A <- matrix(0, nrow = k, ncol = k * p)
    A[1:nrow(Ai[[1]]), 1:ncol(Ai[[1]])] <- Ai[[1]]
    for (i in 2:length(Ai)) {
        if (is.null(Ai[[i]]) == TRUE) {
            Ai[[i]] <- matrix(0, nrow = k, ncol = k)
        }
        else {
            A[1:nrow(Ai[[1]]), 1:ncol(Ai[[1]]) + (i - 1) * ncol(Ai[[1]])] <- Ai[[i]]
        }
    }
    d <- diag((k * p - k))
    d1 <- matrix(0, nrow = nrow(d), ncol = k)
    d <- cbind(d, d1)
    A <- rbind(A, d)
    return(A)
}

# Simulates a VAR process
MultVarSim <- function(k,A1,p,Sigma,n)
  {
Y <- matrix(0,nrow=n+500,ncol=k)
YY <- as.vector(Y)

   for (i in seq(from=(k * p+1),to=(nrow(Y) * k-1),by=k)) {
  u <- as.vector(c(mvrnorm(1,rep(0,k),Sigma),rep(0,k*p-k)))
   YY[(i+k):(i-k*p+1+k)] <- A1%*%YY[(i):(i-k*p+1)]+as.matrix(u,ncol=1)
  
  
   }
YY <- YY[1:(nrow(Y)*k)]
Y <- matrix(YY,ncol=k,byrow=TRUE)
Y <- Y[500:nrow(Y),]
return(Y)
}
# Construct the lag matrix
Zmat2 <- function (Y, p2, k) 
{
    T <- nrow(Y)
    Zcol <- function(q) {
        YPS <- as.vector(t(apply(Y[(1 + q):(p2 + q), ], 2, rev)))
        return(YPS)
    }
    q = 0:(T - 1 - (2 * p2 - p2))
    YPS2 <- sapply(X = q, FUN = Zcol)
    return(list(Z = YPS2, Y = Y[(p2 + 1):nrow(Y), ]))
}

 LambdaGrid<-function (gran1, gran2, Y, Z) 
{
      gamstart = max(t(Y) %*% t(Z))
    gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), 
        length = gran2))
    return(gamm)
}

# Creates beta matrices for each value of lambda
betaCreate <- function(k,p,gran2)
  {
beta=list()
beta1 <- matrix(0,nrow=k,ncol=k*p+1)
for( i in 1:gran2)
  {

    beta[[i]] <- beta1


    }
return(beta)
}

