library(BigVAR)
library(foreach)
library(doParallel)
library(vars)
library(snow)
library(doSNOW)

library(parallel)
cl <- makeCluster(4)

clusterEvalQ(cl,library(BigVAR))
clusterEvalQ(cl,library(vars))
clusterEvalQ(cl,library(MASS))



registerDoSNOW(cl)

YY=list()
for(i in 1:5)
    {
YY[[i]] <-  MultVarSim(k,A1=A,p,.01*diag(k),50)

        }

clusterExport(cl,"YY")
clusterExport(cl,c("p","k"))


system.time(
MSE <- foreach(i=1:5) %dopar%{
## YYY <- MultVarSim(k,A1=AX,p,diag(k),100)
A <- constructModel(YY[[i]],p,"None",c(150,10),RVAR=FALSE,1,cv="Rolling",MN=FALSE,verbose=FALSE)
resH <- cv.BigVAR(A)
}
)   

