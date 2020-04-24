

setClass('sanpled_simset',
         contains='simset',
         representation=list(log.sampling.density='numeric'))

setGeneric('set.up.samples',
           def=function(n.samples,
                        sampling.distribution,
                        method=c('random','lhs','weighted'),
                        prior.simset=NULL,
                        weight.function=NULL,
                        log.prior=NULL,
                        log.likelihood=NULL){
               standardGeneric('set.up.samples')
           })

setMethod('set.up.samples',
          signature(sampling.distribution='Distribution'),
def=function(n.samples,
               sampling.distribution,
               method=c('random','lhs','weighted'),
               prior.simset=NULL,
               weight.function=NULL,
               log.prior=NULL,
               log.likelihood=NULL)
{

})

setGeneric('simulate.samples.from.setup',
           def=function(setup,
                        sample.indices=1:setup@n.sim,
                        cores=parallel::detectCores()){
               standardGeneric('simulate.samples')
           })
setMethod('simulate.samples.from.setup',
          signature(setup='simset_setup'),
def=function(setup,
           sample.indices=1:setup@n.sim,
           cores=parallel::detectCores())
{
    cores = min(cores, parallel::detectCores())

    n.sim = length(sample.indices)
    chunk.sizes = rep(floor(n.sim/cores), length(cores)) + as.numeric((n.sim%%cores) >= (1:cores))
    chunk.index.offsets = c(0,cumsum(chunk.sizes))
    chunk.indices = lapply(1:cores, function(core){
        sample.indices[chunk.indices[core]:chunk.indices[core+1]]
    })

    component.simsets = parallel.lapply(chunk.indices, function(indices){
        one.simset = new('sampling_simset', n.sim=length(chunk.indices),
                         parameter.names = setup@parameter.names,
                         n.parameters = setup@n.parameters)

        one.simset@parameters = setup@parameters[chunk.indices,]
        one.simset@log.sampling.density = setup@weights[chunk.indices]

        one.simset@simulations = lapply(chunk.indices, function(i){
            setup@simulation.function(parameters[i,])
        })
    })

    if (length(component.simsets)==1)
        component.simsets[[1]]
    else
        merge(component.simsets)
})
