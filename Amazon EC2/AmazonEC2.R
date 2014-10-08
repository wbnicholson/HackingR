# Packages Required:
#AWS.tools
#foreach
#doSNOW
library(AWS.tools)
library(doSNOW)
library(foreach)


# need to re-run these functions from AWS.tools so that they are added to the namespace 
sleep.while.pending <- function(reservation.id,sleep.time=2,verbose=TRUE) {
    while(pending.instance(reservation.id)) {
        if(verbose) { cat(".") }
        Sys.sleep(sleep.time)
    }
    if(verbose) { cat("\n") }
}

instances.from.reservation <- function(reservation.id,verbose=FALSE) {
    ec2din(filters=paste("reservation-id",reservation.id,sep="="),verbose=verbose)[[reservation.id]]
}

pending.instance <- function(reservation.id) {
    instances <- instances.from.reservation(reservation.id)

    ## test both state and public dns status
    if(any(is.na(instances[,"instanceState"]) | is.na(instances[,"dnsName"]))) {
        ans <- TRUE
    } else {
        ans <- any(instances[,"instanceState"]=="pending")
    }
    ans
}

# I edited the startCluster function so that it has a "security groups" option.
startCluster <- function (ami, key, instance.count, instance.type,security.groups, verbose = FALSE) 
{
    cmd <- paste("ec2-run-instances", ami, "--show-empty-fields", 
        "--key", key, "--instance-count", instance.count, "--instance-type", 
        instance.type, "--group",security.groups)
    if (verbose) {
        cat("using this cmd:\n")
        print(cmd)
    }
    res <- system(cmd, intern = TRUE)
    reservation <- strsplit(res[[1]], split = "\t")[[1]][-1]
    sleep.while.pending(reservation[1], verbose)
    instances <- instances.from.reservation(reservation[1])
    ans <- list(reservation = reservation, instances = instances)
    class(ans) <- "ec2.cluster"
    ans
}
# start cluster here, ami is based on a "worker node" machine image that I created
# We can only have 20 instances running at once without prior permission from amazon.
cl <- startCluster(ami="ami-xxxx",key="Keyname",instance.count=19,instance.type="t2.micro",security.groups="sg-xxxx",verbose=TRUE)


machines <- cl$instances$dnsName
setDefaultClusterOptions(port=10187)
clust <- makeCluster(machines,type="SOCK")
registerDoSNOW(clust)


##### Do parallel work here

stopCluster(clust)
# stops the SNOW cluster
terminateCluster(cl)
# terminates the EC2 worker instances, so billing stops.  

