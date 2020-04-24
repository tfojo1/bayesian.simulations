#The windows version is adapted largely from:
#https://www.r-bloggers.com/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/
parallel.lapply <- function(X, FUN, cores=parallel::detectCores())
{
    cores = min(cores, length(X), parallel::detectCores())

    if (cores==1)
        lapply(X, FUN)
    else if (.Platform$OS.type=='windows')
    {
        cl <- parallel::makeCluster(cores, outfile="")

        ## Find out the names of the loaded packages
        loaded.package.names <- c(
            ## Base packages
            sessionInfo()$basePkgs,
            ## Additional packages
            names( sessionInfo()$otherPkgs ))

        tryCatch( {

            # Export everything from all ancestor environments
            env.list = list(parent.env(environment()))
            while (!identical(env.list[[1]], globalenv()))
                env.list = c(list(parent.env(env.list[[1]])), env.list)

            for (one.env in env.list)
                parallel::clusterExport(cl,
                                        ls(env=one.env, all.names=TRUE),
                                        envir=one.env)

            # Export just X and FUN from this environment (this function)
            parallel::clusterExport(cl,
                                    c('X', 'FUN'),
                                    envir=environment())

            ## Load the libraries on all the clusters
            ## N.B. length(cl) returns the number of clusters
            parallel::parLapply( cl, 1:length(cl), function(xx){
                lapply(loaded.package.names, require, character.only=T)
            })

            ## Run the lapply in parallel
            return(
                parallel::parLapply(cl=cl, X=X, fun=FUN)
            )

        }, finally = {
            ## Stop the cluster
            parallel::stopCluster(cl)
        })
    }
    else
        parallel::mclapply(X, FUN, mc.cores=cores)
}
