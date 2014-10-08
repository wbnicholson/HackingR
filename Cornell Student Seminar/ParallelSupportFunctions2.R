
# functions to extract relevant information from the s4 object class
returnMSFE <- function(foo)
    {
    MSE <- foo@OOSMSFE   
return(MSE)
}

# Extract MSFE 
returnMeanMSFE <- function(MSE)
    {
return(mean((sapply(lapply(MSE,returnMSFE),"[[",1))))


    }

returnSDMSFE <- function(MSE)
    {
return(sd((sapply(lapply(MSE,returnMSFE),"[[",1))))


    }



returnSD <- function(foo)
    {
se <- foo@seoosmsfe

return(se)
}

# Extract Standard Error
returnSDAgg <- function(se)
    {
        
return(mean((sapply(lapply(se,returnSD),"[[",1))))

        }


returnAIC <- function(foo)
    {
return(foo@AICMSFE)
        }

returnMeanAIC <- function(MSE)
    {
return(mean((sapply(lapply(MSE,returnAIC),"[[",1))))

    }
returnSDAIC <- function(MSE)
    {
return(sd((sapply(lapply(MSE,returnAIC),"[[",1))))

    }


returnAICSD <- function(foo)
    {
return(foo@AICSD)
        }

returnAICSDMean <- function(se)
    {

return(mean((sapply(lapply(se,returnAICSD),"[[",1))))

        }

returnRW <- function(foo)
    {
return(foo@RWMSFE)
        }

returnMeanRW <- function(MSE)
    {
return(mean((sapply(lapply(MSE,returnRW),"[[",1))))

    }
returnSDRW <- function(MSE)
    {
return(sd((sapply(lapply(MSE,returnRW),"[[",1))))

    }




returnRWSD <- function(foo)
    {
return(foo@RWSD)
        }

returnRWSDMean <- function(se)
    {

return(mean((sapply(lapply(se,returnRWSD),"[[",1))))

       }

returnUncondMean <- function(foo)
    {
return(foo@MeanMSFE)
        }

returnMeanMean <- function(MSE)
    {
return(mean((sapply(lapply(MSE,returnUncondMean),"[[",1))))

    }

returnSDMeanMean <- function(MSE)
    {
return(sd((sapply(lapply(MSE,returnUncondMean),"[[",1))))

    }


returnUcondSD <- function(foo)
    {
return(foo@MeanSD)
        }

returnUncondSDMean <- function(se)
    {

return(mean((sapply(lapply(se,returnUcondSD),"[[",1))))

       }

