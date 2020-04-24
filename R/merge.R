
##------------------##
##-- SERIAL MERGE --##
##------------------##

#'@title Merge MCMC objects in serial
#'
#'@description Merges two or more mcmcsim objects by treating the chains in each a serial chains (ie, extensions of the previous chains)
#'
#'@param ... One or more mcmcsim objects or lists containing only mcmcsim objects
#'
#'@export
mcmc.merge.serial <- function(...)
{
    arguments = list(...)
    mcmc.objects = list()

    for (elem in arguments)
    {
        if (is(elem, 'mcmcsim'))
            mcmc.objects = c(mcmc.objects, list(elem))
        else if (is(elem, 'list'))
        {
            not.mcmcsim = !sapply(elem, function(sub.elem){is(sub.elem, 'mcmcsim')})
            if (any(not.mcmcsim))
                stop("The arguments passed in ... must be either mcmcsim objects or lists of mcmcsim objects")
            mcmc.objects = c(mcmc.objects, elem)
        }
        else
            stop("The arguments passed in ... must be either mcmcsim objects or lists of mcmcsim objects")
    }

    # If length 1, just return
    if (length(mcmc.objects)==1)
        return(mcmc.objects[[1]])

    # If length 0, throw an error
    if (length(mcmc.objects)==0)
        stop("mcmc.objects is an empty list")

    # Check chain length
    chains.per = sapply(mcmc.objects, function(mcmc){mcmc@n.chains})
    if (length(unique(chains.per))>1)
        stop("All elements of mcmc.objects must have the same number of chains")

    # Check that controls are mergable
    c1 = mcmc.objects[[1]]@control
    sapply(2:length(mcmc.objects), function(i){
        errors = check.merge.controls(c1, mcmc.objects[[i]]@control)
        if (length(errors)>0)
            stop(paste0("Unable to merge the ",
                        get.ordinal(i),
                        " mcmc object with those prior. The controls are incompatible: ",
                        paste0(errors, collapse=', ')))
    })

    # MCMC variables
    n.iter = sum(sapply(mcmc.objects, function(mcmc){mcmc@n.iter}))
    n.chains = mcmc.objects[[1]]@n.chains
    n.var = c1@n.var

    # Set up RV data structure
    rv = create.skeleton.mcmc(c1,
                              n.iter=n.iter,
                              n.chains=n.chains)

    #Loop through and merge
    sim.counter = 0
    iter.inc = 0

    rv@n.accepted.in.burn[] = as.numeric(0)
    rv@run.times.in.burn[] = as.numeric(0)
    rv@total.run.time[] = as.numeric(0)

    for (mcmc in mcmc.objects)
    {
        # first, check for repeated simulations
        first.sim.is.repeat.by.chain = iter.inc > 0 &
            sapply(1:n.chains, function(chain){
                if (iter.inc==1)
                    all(mcmc@samples[chain,1,] == rv@samples[chain,1,])
                else if (n.var==1)
                    any(rv@samples[chain,1:iter.inc,1] == mcmc@samples[chain,1,1])
                else if (iter.inc>0)
                    any(apply(rv@samples[chain,1:iter.inc,], 'iteration', function(ss){all(ss==mcmc@samples[chain,1,])}))
                else
                    F
            })

        #Simulations
        sims.to.remove = mcmc@simulation.indices[first.sim.is.repeat.by.chain,1]
        rv@simulations = c(rv@simulations, mcmc@simulations[setdiff(1:length(mcmc@simulations), sims.to.remove)])

        #Sim indices
        rv@simulation.indices[,iter.inc + 1:mcmc@n.iter] = mcmc@simulation.indices + sim.counter - as.numeric(first.sim.is.repeat.by.chain)
        sim.counter = length(rv@simulations)

        #Samples
        rv@samples[,iter.inc + 1:mcmc@n.iter,] = mcmc@samples

        #Likelihood and prior
        rv@log.likelihoods[,iter.inc + 1:mcmc@n.iter] = mcmc@log.likelihoods
        rv@log.priors[,iter.inc + 1:mcmc@n.iter] = mcmc@log.priors

        #n.accepted
        rv@n.accepted[,iter.inc + 1:mcmc@n.iter,] = mcmc@n.accepted

        #run.times
        rv@run.times[,iter.inc + 1:mcmc@n.iter] = mcmc@run.times

        #first step in iter
        rv@first.step.for.iter[,iter.inc + 1:mcmc@n.iter] = mcmc@first.step.for.iter

        #update iter.inc
        iter.inc = iter.inc + mcmc@n.iter

        #add to the cumulative tally states
        rv@n.accepted.in.burn = rv@n.accepted.in.burn + mcmc@n.accepted.in.burn
        rv@run.times.in.burn = rv@run.times.in.burn + mcmc@run.times.in.burn
        rv@total.run.time = rv@total.run.time + mcmc@total.run.time
    }

    #Add the chain state to the rv
    rv@chain.states = mcmc.objects[[length(mcmc.objects)]]@chain.states

    # Return
    rv
}

##--------------------##
##-- PARALLEL MERGE --##
##--------------------##


#'@title Merge MCMC objects in parallel
#'
#'@description Merges two or more mcmcsim objects by treating the chains in each as parallel chains
#'
#'@export
mcmc.merge.parallel <- function(...)
{
    arguments = list(...)
    mcmc.objects = list()

    for (elem in arguments)
    {
        if (is(elem, 'mcmcsim'))
            mcmc.objects = c(mcmc.objects, list(elem))
        else if (is(elem, 'list'))
        {
            not.mcmcsim = !sapply(elem, function(sub.elem){is(sub.elem, 'mcmcsim')})
            if (any(not.mcmcsim))
                stop("The arguments passed in ... must be either mcmcsim objects or lists of mcmcsim objects")
            mcmc.objects = c(mcmc.objects, elem)
        }
        else
            stop("The arguments passed in ... must be either mcmcsim objects or lists of mcmcsim objects")
    }

    # If length 1, just return
    if (length(mcmc.objects)==1)
        return(mcmc.objects[[1]])

    # If length 0, throw an error
    if (length(mcmc.objects)==0)
        stop("mcmc.objects is an empty list")

    # Check n.iter
    iter.per = sapply(mcmc.objects, function(mcmc){mcmc@n.iter})
    if (length(unique(iter.per))>1)
        stop("All elements of mcmc.objects must have the same number of iterations")

    # Check that controls are mergable
    c1 = mcmc.objects[[1]]@control
    sapply(2:length(mcmc.objects), function(i){
        errors = check.merge.controls(c1, mcmc.objects[[i]]@control)
        if (length(errors)>0)
            stop(paste0("Unable to merge the ",
                        get.ordinal(i),
                        " mcmc object with those prior. The controls are incompatible: ",
                        paste0(errors, collapse=', ')))
    })

    # MCMC variables
    n.iter = mcmc.objects[[1]]@n.iter
    n.chains = sum(sapply(mcmc.objects, function(mcmc){mcmc@n.chains}))
    n.var = c1@n.var

    # Set up RV data structure
    rv = create.skeleton.mcmc(c1,
                              n.iter=n.iter,
                              n.chains=n.chains)

    #Loop through and merge
    sim.counter = 0
    chain.counter = 0

    rv@n.accepted.in.burn[] = as.numeric(0)
    rv@run.times.in.burn[] = as.numeric(0)
    rv@total.run.time[] = as.numeric(0)

    for (mcmc in mcmc.objects)
    {
        chains = chain.counter + 1:mcmc@n.chains
        chain.counter = chain.counter + mcmc@n.chains

        #No need to check for repeated simulations with parallel merge
        #Simulations
        rv@simulations = c(rv@simulations, mcmc@simulations)

        #Sim indices
        rv@simulation.indices[chains,] = mcmc@simulation.indices[1:mcmc@n.chains,] + sim.counter
        sim.counter = length(rv@simulations)

        #Samples
        rv@samples[chains,,] = mcmc@samples[1:mcmc@n.chains,,]

        #Likelihood and prior
        rv@log.likelihoods[chains,] = mcmc@log.likelihoods[1:mcmc@n.chains,]
        rv@log.priors[chains,] = mcmc@log.priors[1:mcmc@n.chains,]

        #n.accepted
        rv@n.accepted[chains,,] = mcmc@n.accepted[1:mcmc@n.chains,,]

        #run.times
        rv@run.times[chains,] = mcmc@run.times[1:mcmc@n.chains,]

        #first step in iter
        rv@first.step.for.iter[chains,] = mcmc@first.step.for.iter[1:mcmc@n.chains,]


        #add to the cumulative tally states
        rv@n.accepted.in.burn[chains,] = mcmc@n.accepted.in.burn[1:mcmc@n.chains,]
        rv@run.times.in.burn[chains] = mcmc@run.times.in.burn
        rv@total.run.time[chains] = mcmc@total.run.time

        #chain states
        rv@chain.states[chains] = mcmc@chain.states
    }

    # Return
    rv
}

##-----------------------------##
##-- CHECKING IF OK TO MERGE --##
##-----------------------------##

do.check.merge.controls <- function(c1, c2, for.serial.merge)
{
    errors = character()
    # Class
    if (!all(class(c1)==class(c2)))
        errors = c(errors, "Control classes do not match")

    # Vars
    if (c1@n.var != c2@n.var || any(c1@var.names != c2@var.names))
        errors = c(errors, "Variable names do not match")

    # Method
    if (any(c1@method != c2@method))
        errors = c(errors, "Methods are not the same")

    # Thin
    if (c1@thin != c2@thin)
        errors = c(errors, "Thinning is not the same")

    # Burn
    if (c1@burn != c2@burn)
        errors = c(errors, "Burn is not the same")

    # Sample steps
    if (length(c1@sample.steps) != length(c2@sample.steps))
        errors = c(errors, "Different number of sample steps")
    else if (any(c1@sample.steps != c2@sample.steps))
        errors = c(errors, "Sample steps are not the same")

    # Return
    errors
}
