
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


setClass('mcmcsim_control',
         representation=list(var.names='character',
                             n.var='integer',
                             method='character',
                             thin='integer',
                             burn='integer',
                             simulation.function='function',
                             pass.chain.to.simulation.function='logical',
                             pass.iteration.to.simulation.function='logical',
                             log.prior.distribution='function',
                             log.likelihood='function',
                             transformations='list',
                             lower.bounds='numeric',
                             upper.bounds='numeric',
                             sample.steps='character'))

setClass('mcmcsim',
         representation=list(control='mcmcsim_control',
                             var.names='character',
                             n.var='integer',
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

##-------------------------------------------------------------------##
##-- CLASS METHODS THAT *MAY* BE OVERRIDDEN AT THE SUB-CLASS LEVEL --##
##-------------------------------------------------------------------##

#'@title Check whether mcmc objects resulting from two mcmcsim_control objects can be merged
#'@return If the two controls are merge-able, returns an empty character vector. Otherwise, returns a vector of the reasons why the two cannot be merged
#'@export
setGeneric("check.merge.controls", function(c1, c2, for.serial.merge=T)
{
    standardGeneric('check.merge.controls')
})
setMethod("check.merge.controls",
          signature(c1='mcmcsim_control', c2='mcmcsim_control'),
          def = function(c1, c2, for.serial.merge=T)
{
    do.check.merge.controls(c1, c2, for.serial.merge)
})



##--------------------------##
##-- PRIVATE CONSTRUCTORS --##
##--------------------------##


# Creates a skeleton MCMC object in which to put samples
create.skeleton.mcmc <- function(control,
                                 n.iter,
                                 n.chains=1,
                                 class.name='mcmcsim',
                                 ...)
{
    # Set array names
    if (n.iter==0)
        iterations.name = character()
    else
        iterations.name = 1:n.iter

    sample.dimnames = list(chain=1:n.chains, iteration=iterations.name, variable=control@var.names)
    sim.index.dimnames = sample.dimnames[1:2]
    n.accepted.dimnames = list(chain=1:n.chains, iteration=iterations.name, step=control@sample.steps)

    # Make data structures
    simulation.indices = matrix(as.integer(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    samples = array(as.numeric(NA), dim=sapply(sample.dimnames, length), dimnames=sample.dimnames)
    log.likelihoods = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    log.priors = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    n.accepted = array(as.integer(NA), dim=sapply(n.accepted.dimnames, length), dimnames=n.accepted.dimnames)
    run.times = matrix(as.numeric(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)
    first.step.for.iter = matrix(as.integer(NA), nrow=n.chains, ncol=n.iter, dimnames=sim.index.dimnames)

    run.times.in.burn = rep(as.numeric(0), n.chains)
    n.accepted.in.burn = matrix(as.integer(0), nrow=n.chains, ncol=length(control@sample.steps),
                                dimnames=list(chain=1:n.chains, step=control@sample.steps))
    total.run.time = rep(as.numeric(NA), n.chains)
    simulations = list()

    # Package and return
    new(Class = class.name,
        control=control,
        var.names=control@var.names,
        n.var=control@n.var,
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
}

create.mcmcsim.control <- function(class.name,
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

    if (thin < 1)
        stop("thin must be >= 1")

    if (burn<0)
        stop("burn must be >= 0")

    new(Class=class.name,
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
}


create.new.mcmcsim.chain.state <- function(class.name, current.parameters, ...)
{
    new(Class=class.name,
        current.parameters=current.parameters,
        n.unthinned.after.burn=as.integer(0),
        n.unthinned.burned=as.integer(0),
        ...)
}
