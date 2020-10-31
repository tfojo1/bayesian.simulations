
# In general, use of MCMC here is going to go as follows
# Create an MCMC control object
# Run one or more chains
# Merge the chains into a single MCMC object
# Extend chains if desired
# Extract the posterior distribution
# Profit


##------------------------------##
##-- MASTER CLASS DEFINITIONS --##
##------------------------------##

#'@importClassesFrom distributions Distribution
setClassUnion('Distribution_or_function',
              members=c('Distribution',
                        'function'))

setClass('mcmcsim_control',
         representation=list(var.names='character',
                             n.var='integer',
                             method='character',
                             thin='integer',
                             burn='integer',
                             simulation.function='function',
                             pass.chain.to.simulation.function='logical',
                             pass.iteration.to.simulation.function='logical',
                             log.prior.distribution='Distribution_or_function',
                             log.likelihood='function',
                             transformations='list',
                             lower.bounds='numeric',
                             upper.bounds='numeric',
                             sample.steps='character'))

setClass('mcmcsim',
         representation=list(var.names='character',
                             n.var='integer',
                             method='character',
                             thin='integer',
                             burn='integer',
                             transformations='list',
                             lower.bounds='numeric',
                             upper.bounds='numeric',
                             sample.steps='character',

                             simulations='list',
                             simulation.indices='matrix',
                             samples='array',
                             log.likelihoods='matrix',
                             log.priors='matrix',
                             n.chains='integer',
                             chain.states='list',
                             n.iter='integer',
                             n.accepted='array',
                             run.times='matrix', #indexed [chain,iter]
                             first.step.for.iter='matrix',
                             n.accepted.in.burn='matrix',
                             run.times.in.burn='numeric', #indexed [chain]
                             total.run.time='numeric'))

setClass('mcmcsim_chainstate',
         representation=list(
             current.parameters='numeric',
             n.unthinned.after.burn='integer',
             n.unthinned.burned='integer',
             n.accepted='integer',
             run.time='numeric',
             first.step.for.iter='integer'
         ))


##-------------------------------------------------------------------##
##-- CLASS METHODS THAT WILL BE IMPLEMENTED AT THE SUB-CLASS LEVEL --##
##--                                                               --##
##--        (To be called internally by the MCMC package -         --##
##--                    not directly by users)                     --##
##--                                                               --##
##-------------------------------------------------------------------##

setGeneric("run.single.chain", function(control,
                                        chain.state,
                                        n.iter,
                                        update.frequency=floor(n.iter/10),
                                        update.detail='med',
                                        initial.sim=NULL,
                                        total.n.iter=n.iter,
                                        prior.n.iter=0,
                                        prior.n.accepted=0,
                                        prior.run.time=0,
                                        return.current.sim=F,
                                        chain,
                                        output.stream)
{
    standardGeneric("run.single.chain")
})

setGeneric("create.initial.chain.state", function(control,
                                                  start.parameters)
{
    standardGeneric("create.initial.chain.state")
})



##--------------------------##
##-- PRIVATE CONSTRUCTORS --##
##--------------------------##


setMethod("initialize",
          signature(.Object='mcmcsim'),
          def=function(.Object,
                       n.iter,
                       n.chains,
                       var.names=model@var.names,
                       n.var=model@n.var,
                       method=model@method,
                       thin=model@thin,
                       burn=model@burn,
                       transformations=model@transformations,
                       lower.bounds=model@lower.bounds,
                       upper.bounds=model@upper.bounds,
                       sample.steps=model@sample.steps,
                       model=NULL)
{
    # Set array names
    if (n.iter==0)
        iterations.name = character()
    else
        iterations.name = (1:n.iter) * thin + burn

    sample.dimnames = list(chain=1:n.chains, iteration=iterations.name, variable=var.names)
    sim.index.dimnames = sample.dimnames[1:2]
    n.accepted.dimnames = list(chain=1:n.chains, iteration=iterations.name, step=sample.steps)

    # Make data structures
    simulation.indices = matrix(as.integer(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    samples = array(as.numeric(NA), dim=sapply(sample.dimnames, length), dimnames=sample.dimnames)
    log.likelihoods = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    log.priors = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    n.accepted = array(as.integer(NA), dim=sapply(n.accepted.dimnames, length), dimnames=n.accepted.dimnames)
    run.times = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    first.step.for.iter = matrix(as.integer(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)

    run.times.in.burn = rep(as.numeric(0), n.chains)
    n.accepted.in.burn = matrix(as.integer(0), nrow=n.chains, ncol=length(sample.steps),
                                dimnames=list(chain=1:n.chains, step=sample.steps))
    total.run.time = rep(as.numeric(NA), n.chains)
    simulations = list()

    # Package and return
    callNextMethod(.Object,
                   var.names=var.names,
                   n.var=n.var,
                   method=method,
                   thin=thin,
                   burn=burn,
                   transformations=transformations,
                   lower.bounds=lower.bounds,
                   upper.bounds=upper.bounds,
                   sample.steps=sample.steps,

                   simulations=simulations,
                   simulation.indices=simulation.indices,
                   samples=samples,
                   log.likelihoods=log.likelihoods,
                   log.priors=log.priors,
                   n.chains=as.integer(n.chains),
                   n.iter=as.integer(n.iter),
                   n.accepted=n.accepted,
                   first.step.for.iter=first.step.for.iter,
                   run.times=run.times,
                   n.accepted.in.burn=n.accepted.in.burn,
                   run.times.in.burn=run.times.in.burn,
                   total.run.time=total.run.time)
})

setMethod("initialize",
          signature(.Object='mcmcsim_control'),
          def=function(.Object,
                       var.names,
                       method,
                       simulation.function,
                       pass.chain.to.simulation.function,
                       pass.iteration.to.simulation.function,
                       log.prior.distribution,
                       log.likelihood,
                       thin,
                       burn,
                       transformations,
                       sample.steps,
                       ...)
{
    #process transformations
    transformations = process.transformations(var.names, transformations)

    #do we need to check sample steps against var names?

    if (thin < 1)
        stop("thin must be >= 1")

    if (burn<0)
        stop("burn must be >= 0")

    callNextMethod(.Object,
                   var.names=var.names,
                   n.var=length(var.names),
                   method=method,
                   simulation.function=simulation.function,
                   pass.chain.to.simulation.function=pass.chain.to.simulation.function,
                   pass.iteration.to.simulation.function=pass.iteration.to.simulation.function,
                   log.prior.distribution=log.prior.distribution,
                   log.likelihood=log.likelihood,
                   thin=thin,
                   burn=burn,
                   transformations=transformations,
                   sample.steps=sample.steps,
                   ...
    )
})

setMethod("initialize",
          signature(.Object='mcmcsim_chainstate'),
          def=function(.Object,
                       current.parameters,
                       ...)
{
    callNextMethod(.Object,
                   current.parameters=current.parameters,
                   n.unthinned.after.burn=as.integer(0),
                   n.unthinned.burned=as.integer(0),
                   ...)
})



##-----------------##
##-- SUBSET MCMC --##
##-----------------##

#'@title Subset an MCMC object
#'
#'@description Subset an MCMC object by either burning additional iterations of the front, further thinning iterations, selecting a subset of chains, or some combination of the three
#'
#'@param .Object The mcmc.sim object
#'@param chains The chains to subset
#'@param additional.burn Additional number of iterations to burn off the front (this represents the number of iterations AFTER the thinning which has already been done to the mcmc.sim object but BEFORE applying additional.thin)
#'@param additional.thin The factor by which to further thin iterations
#'
#'@return An mcmc.sim object
#'
#'@export
setGeneric("mcmc.subset",
           def=function(.Object,
                        chains=1:.Object@n.chains,
                        additional.burn=0,
                        additional.thin=1){
               standardGeneric("mcmc.subset")
           })
setMethod("mcmc.subset",
          signature(.Object='mcmcsim'),
def=function(.Object,
           chains=1:.Object@n.chains,
           additional.burn=0,
           additional.thin=1)
{
    ## UPDATE THE MCMC OBJECT
    # Save some initial values
    orig.thin = .Object@thin
    orig.burn = .Object@burn
    orig.n.iter = .Object@n.iter
    orig.n.chains = .Object@n.chains

    # Call the sub-function and set up indices
    subsetted = do.subset.mcmc(mcmc=.Object,
                               chains=chains,
                               additional.burn=additional.burn,
                               additional.thin=additional.thin)
    new.n.iter = subsetted$n.iter
    new.n.chains = length(chains)
    n.newly.burned = as.integer(additional.burn * orig.thin)
    new.burn = .Object@burn + n.newly.burned
    new.thin = as.integer(.Object@thin * additional.thin)

    keep.indices = (1:orig.n.iter)[subsetted$keep.mask]
    first.in.thin.group.indices = (1:orig.n.iter)[subsetted$first.in.thin.group.mask]
    burned.mask = (1:orig.n.iter)<=additional.burn
    indices.after.last.keep = (1:orig.n.iter)[(1:orig.n.iter)>max(keep.indices)]

    collapsing.matrix = sapply(1:new.n.iter, function(i){
        col = rep(0, orig.n.iter)
        col[first.in.thin.group.indices[i]:keep.indices[i]] = 1
        col
    })

    # Set up new dim.names
    iteration.names = new.burn + (1:new.n.iter) * new.thin
    dim.names.2d = list(chain=1:new.n.chains, iteration=iteration.names)
    dim.names.3d = list(chain=1:new.n.chains, iteration=iteration.names, variable=.Object@var.names)

    # Update samples
    .Object@samples = subsetted$samples

    # Update n.iter and n.chains
    .Object@n.iter = new.n.iter
    .Object@n.chains = new.n.chains

    # Update thin and burn
    .Object@burn = new.burn
    .Object@thin = new.thin

    # Update simulations and simulation indices
    new.to.old.sim.indices = unique(as.numeric(t(subsetted$simulation.indices)))
    n.new.sims = length(new.to.old.sim.indices)
    .Object@simulations = .Object@simulations[new.to.old.sim.indices]
    .Object@simulation.indices = subsetted$simulation.indices
    .Object@simulation.indices[] = sapply(.Object@simulation.indices[], function(old.i){
        (1:n.new.sims)[new.to.old.sim.indices==old.i]
    })

    # Update log.likelihoods and log.priors
    .Object@log.likelihoods = .Object@log.likelihoods[chains, keep.indices]
    .Object@log.priors = .Object@log.priors[chains, keep.indices]
    dim(.Object@log.likelihoods) = dim(.Object@log.priors) = sapply(dim.names.2d, length)
    dimnames(.Object@log.likelihoods) = dimnames(.Object@log.priors) = dim.names.2d

    # Update n.accepted.in.burn
    if (additional.burn>0)
    {
        newly.burned = .Object@n.accepted[chains,keep.indices,]
        if (sum(burned.mask)==1 && new.n.chains==1)
            newly.burned = sum(newly.burned)
        else if (new.n.chains==1)
            newly.burned = colSums(newly.burned)
        else if (sum(burned.mask)>1)
            newly.burned = apply(newly.burned, c('chain','step'), sum)

        .Object@n.accepted.in.burn = .Object@n.accepted.in.burn + newly.burned
    }

    # Update n.accepted
    #  save this for updating the chain states
    accepted.after.last.iter = sapply(.Object@sample.steps, function(step){
        sapply(1:orig.n.chains, function(chain){
            sum(.Object@n.accepted[chain,indices.after.last.keep,step])
        })
    })
    dim(accepted.after.last.iter) = c(chain=orig.n.chains, step=length(.Object@sample.steps))

    dim.names.accepted = list(chain=1:new.n.chains, iteration=iteration.names, step=.Object@sample.steps)

    .Object@n.accepted = sapply(.Object@sample.steps, function(step){
        .Object@n.accepted[,,step] %*% collapsing.matrix
    })
    dim(.Object@n.accepted) = sapply(dim.names.accepted, length)
    dimnames(.Object@n.accepted) = dim.names.accepted

    # Update run.times.in.burn
    run.time.after.last.iter = sapply(1:orig.n.chains, function(chain){
        sum(.Object@run.times[chain, indices.after.last.keep])
    })

    if (additional.burn>0)
    {
        if (sum(burned.mask)==1)
            newly.burned = .Object@run.times[chains, burned.mask]
        else if (new.n.chains==1)
            newly.burned = sum(.Object@run.times[chains, burned.mask])
        else
            newly.burned = rowSums(.Object@run.times[chains, burned.mask])

        .Object@run.times.in.burn = .Object@run.times.in.burn[chains] + newly.burned
    }

    # Update run.times
    .Object@run.times = .Object@run.times %*% collapsing.matrix
    dim(.Object@run.times) = sapply(dim.names.2d, length)
    dimnames(.Object@run.times) = dim.names.2d

    # Update first.step.for.iter
    .Object@first.step.for.iter = .Object@first.step.for.iter[chains, first.in.thin.group.indices]
    dim(.Object@first.step.for.iter) = sapply(dim.names.2d, length)
    dimnames(.Object@first.step.for.iter) = dim.names.2d

    # Update total.run.time
    .Object@total.run.time = .Object@total.run.time[chains]

    ## UPDATE THE CHAIN STATES
    .Object@chain.states = lapply(chains, function(chain){
        chain.state = .Object@chain.states[[chain]]

        chain.state@n.unthinned.after.burn = chain.state@n.unthinned.after.burn - n.newly.burned
        chain.state@n.unthinned.burned = chain.state@n.unthinned.burned + n.newly.burned

        chain.state@n.accepted = chain.state@n.accepted + accepted.after.last.iter[chain,]
        chain.state@run.time = chain.state@run.time + run.time.after.last.iter[chain]

        chain.state
    })

    .Object
})
