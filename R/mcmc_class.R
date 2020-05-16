
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
                             run.times='matrix',
                             first.step.for.iter='matrix',
                             n.accepted.in.burn='matrix',
                             run.times.in.burn='numeric',
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
                       control,
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
        iterations.name = 1:n.iter

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
