#Pre-parallelization Benchmarking
# require package microbenchmark, ggplot 2 to recreate the plot

source("LassoVARAlgorithms.R")
source("LassoVARSupportFunctions.R")

# microbenchmark is used to show performance differences
## install.packages('microbenchmark')
library(microbenchmark)
library(MASS)

p=2;k=4

# Start by generating some data
Ai=list()
Ai[[1]] <- matrix(runif(16,-1,1)*rbinom(64,1,.5),nrow=4,ncol=4)
Ai[[2]] <- matrix(0,nrow=4,ncol=4)
A <- VarptoVar1(Ai,p,k)
# check to ensure stationarity
while(max(Mod(eigen(A)$values))>1)
{
Ai=list()
Ai[[1]] <- matrix(runif(16,-1,1)*rbinom(64,1,.5),nrow=4,ncol=4)
Ai[[2]] <- matrix(0,nrow=4,ncol=4)
A <- VarptoVar1(Ai,p,k)
}

nsims=100

# simulate 100 observations then create lag matrices
YYY <- MultVarSim(k,A,p,diag(k),nsims)

ZZZ <- Zmat2(YYY,p,k)
Z <- ZZZ$Z
Y <- ZZZ$Y


gran2=50 # creates 50 penalty parameters
#generate penalty parameters
gamm <- LambdaGrid(40,gran2,Y,Z)
# generate coefficient matrices
# lists for non-optimized methods
B <- betaCreate(k,p=2,gran2)

# array for parallel method
B1a <- array(0,dim=c(k,k*p+1,gran2))
# microbenchmarking

op <- microbenchmark(
B1 <- lassoVAR(B,Z,Y,gamm,1e-3)
  ,
  B2 <- lassoVARVec(B,Z,Y,gamm,1e-3)
  ,
B3 <- lassoVARCPP(B,Z,Y,gamm,1e-3)
,
B4 <- lassoVARCPPPar(B1a,Z,Y,gamm,1e-3)
    
  ,times=50L)


plot(op,main="Computational Time over 50 iterations",xlab="Procedure",names=c(" Three For loops","Vectorized","Exported to Rcpp"))



