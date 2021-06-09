
##----------------------------------------##
##-- CLASS DEFINITIONS and CONSTRUCTORS --##
##----------------------------------------##

setClassUnion('character_or_NULL', members=c('NULL','character'))

setClass('parameter_set',
         representation=list(parameter.names='character_or_NULL',
                             n.parameters='integer',
                             parameters='matrix',
                             n.sim='integer',
                             weights='numeric',
                             parameter.transformations='list',
                             parameter.lower.bounds='numeric',
                             parameter.upper.bounds='numeric'
         ))

setMethod('initialize',
          signature(.Object='parameter_set'),
def=function(.Object,
           n.sim=NULL,
           parameter.names=NULL,
           n.parameters=if(is.null(parameter.names)) NA else length(parameter.names),
           parameter.values=NULL,
           weights=1,
           parameter.transformations=NULL)
{
    # Initial check parameter names
    if (!is.null(parameter.names))
        parameter.names = as.character((parameter.names))
    if (length(parameter.names)==0)
        parameter.names = NULL


    if (!is.null(parameter.values))
    {
        if (!is(parameter.values, 'matrix') && !is(parameter.values, 'array'))
            stop("parameter.values must be a matrix or 2-dimensional array")
        if (length(dim(parameter.values)) != 2)
            stop("parameter.values must be a matrix or 2-dimensional array")
        if (!is(parameter.values[1,1], 'numeric') && !is(parameter.values[1,1], 'integer'))
            stop('parameter.values must contain numeric values only')

        n.sim = dim(parameter.values)[1]
        n.parameters = dim(parameter.values)[2]

        if (is.null(parameter.names))
            parameter.names = dimnames(parameter.values)[[2]]
    }
    else if (!is.null(parameter.names))
        n.parameters = length(parameter.names)
    else if ((!is(n.parameters, 'integer') &&
                !(is(n.parameters, 'numeric') && round(n.parameters)==n.parameters)) ||
               length(n.parameters) != 1 ||
               is.na(n.parameters) || n.parameters<1)
        stop("'n.parameters' must be a positive, non-NA integer")


    if ((!is(n.sim, 'integer') &&
         !(is(n.sim, 'numeric') && round(n.sim)==n.sim)) ||
        length(n.sim) != 1 ||
        is.na(n.sim) || n.sim<1)
        stop("'n.sim' must be a positive, non-NA integer")

    if (!is.null(parameter.names))
    {
        if (length(parameter.names) != n.parameters)
            stop('parameter.names must have one value for each column in parameter.values')

        if (any(is.na(parameter.names)))
            stop('parameter.names cannot have NA values')

        if (max(table(parameter.names))>1)
            stop('parameter.names must be unique')

        names(parameter.names) = NULL
    }

    if (is.null(parameter.values))
        parameter.values = matrix(as.numeric(NA), nrow=n.sim, ncol=n.parameters)

    dimnames(parameter.values) = list(NULL, parameter=parameter.names)

    #check weights
    if (!is.null(weights) && !is(weights, 'numeric') && !is(weights, 'integer'))
        stop("weights must must numeric or integer values")
    if (is.null(weights) || length(weights)==0)
        weights = 1
    if (length(weights)==1)
        weights = rep(weights, n.sim)
    if (length(weights) != n.sim)
        stop(paste0("weights must be either a vector of length n.sim (", n.sim,
                    ") or a scalar value"))

    .Object@parameter.names = parameter.names
    .Object@n.parameters = as.integer(n.parameters)
    .Object@parameters = parameter.values
    .Object@n.sim = as.integer(n.sim)
    .Object@weights = as.numeric(weights)

    .Object@parameter.transformations = process.transformations(.Object@parameter.names,
                                                                parameter.transformations,
                                                                var.names.name='parameter.names')
    #for now
#    .Object@parameter.transformations = list()
#    parameter.lower.bounds=
#    parameter.upper.bounds='numeric'

    .Object
})

setMethod('show',
          signature(object='parameter_set'),
          def=function(object){
              if (all(is.na(object@parameters)))
                  stop(paste0("An empty parameter_set with planned ", object@n.sim, " samples of ", object@n.parameters, " parameters"))
              else
                  stop(paste0("A parameter_set with ", object@n.sim, " samples of ", object@n.parameters, " parameters"))
          })

#'@title The simset class
#'
#'@description A class representing a weighted set of simulations
#'
#'@slot simulations A list of simulations
#'@slot weights A numeric weight for each simulation
#'@slot n.sim The number of unique simulations
#'@slot parameters A matrix of parameters for each simulation, where row i corresponds to the parameters used to run the ith simulation, and each column represents a parameter value
#'@slot parameter.names,n.parameters The names and number of the parameters used to run simulations

#'@name simset
#'@rdname simset
#'@aliases simset-class
#'@exportClass simset
#'@export
setClass('simset',
         contains='parameter_set',
         representation=list(
             simulations='list',
             weights='numeric')
         )

setMethod('initialize',
          signature(.Object='simset'),
def=function(.Object,
             n.sim=NULL,
             parameter.names=NULL,
             n.parameters=if(is.null(parameter.names)) NA else length(parameter.names),
             parameter.values=NULL,
             weights=1,
             simulations=NULL,
             parameter.transformations=NULL)
{
    .Object = callNextMethod(.Object,
                             n.sim=n.sim,
                             parameter.names=parameter.names,
                             n.parameters=n.parameters,
                             parameter.values=parameter.values,
                             weights=weights,
                             parameter.transformations = parameter.transformations)

    if (is.null(simulations))
        simulations = lapply(1:n.sim, function(i){NULL})

    if (!is(simulations, 'list') || length(simulations) != .Object@n.sim)
    {
        if (!is.null(parameter.values))
            stop("simulations must be the same length as the number of rows in parameter.values")
        else
            stop(paste0("simulations must be of length n.sim (", n.sim, ')'))
    }

    .Object@simulations = simulations

    .Object
})

setMethod('show',
          signature(object='simset'),
          def=function(object){
              if (all(is.na(object@parameters)) && all(sapply(object@simulations, is.null)))
                  stop(paste0("An empty simset with planned ", object@n.sim, " simulations from ", object@n.parameters, " parameters"))
              else if (all(is.na(object@parameters)))
                  stop(paste0("A simset with planned ", object@n.sim, " simulations from ", object@n.parameters, " parameters"))
              else
                  stop(paste0("A simset with ", object@n.sim, " simulations from ", object@n.parameters, " parameters"))
          })

#'@title Extracts a simset from an object containing simulations
#'
#'@param object The object from which to extract the simset
#'
#'@export
setGeneric('extract.simset', function(object, ...){
    standardGeneric('extract.simset')
})
setMethod('extract.simset',
          signature(object='mcmcsim'),
          def=function(object,
                       chains=1:object@n.chains,
                       additional.burn=0,
                       additional.thin=1)
          {
              #Check arguments
              chains = check.chains(mcmc=object, chains=chains)

              if (additional.burn < 0 || additional.burn >= object@n.iter)
                  stop("'additional.burn' must be >= 0 and < the number of iteration in the MCMC object (",
                       object@n.iter, ")")
              if (additional.thin < 0 || additional.thin > object@n.iter)
                  stop("'additional.thin' must be >= 0 and <= the number of iteration in the MCMC object (",
                       object@n.iter, ")")


              #Decide what we're going to keep
              keep.mask = ((1:object@n.iter) > additional.burn) &
                  ((1:object@n.iter) %% additional.thin == 0)

              if (sum(keep.mask)==0)
                  stop(paste0("No simulations are available after applying an additional burn of ",
                              additional.burn, " and an additional thin of ", additional.thin))

              dim.names = list(chain=chains, iteration=(1:object@n.iter)[keep.mask], variable=object@var.names)

              raw.parameters = object@samples[chains, keep.mask, ]
              dim(raw.parameters) = sapply(dim.names, length)
              dimnames(raw.parameters) = dim.names

              raw.simulation.indices = object@simulation.indices[chains, keep.mask]
              dim(raw.simulation.indices) = sapply(dim.names[1:2], length)
              dimnames(raw.simulation.indices) = dim.names[1:2]

              keep.simulation.indices = unique(as.integer(raw.simulation.indices))

              new('simset',
                  parameter.names=object@var.names,
                  simulations=object@simulations[keep.simulation.indices],
                  parameter.values = t(sapply(keep.simulation.indices, function(i){
                      mask = raw.simulation.indices==i
                      params = raw.parameters[mask]
                      dim(params) = c(sum(mask), object@n.var)
                      params[1,]
                  })),
                  weights=sapply(keep.simulation.indices, function(i){
                      sum(raw.simulation.indices==i)
                  }),
                  n.sim=length(keep.simulation.indices),
                  n.parameters=object@n.var,
                  parameter.transformations=object@transformations)
               #   parameter.lower.bounds=object@lower.bounds,
               #   parameter.upper.bounds=object@upper.bounds)
          })


#'@title Extract a distribution of simulation results
#'
#'@description Gets a probability distribution over the results of a function applied to each simulation or its parameters or both
#'
#'@param simset A \link{simset} object
#'@param fn A function to be applied. If pass.to.fn=='simulation', the function's first parameters should be 'sim' (a simulation objet). If pass.to.fn=='parameters', the function's first argument should be 'parameters' (a vector of parameters). If pass.to.fn=='both', the function's first two parameters should be 'sim' and 'parameters'
#'@param pass.to.fn What to pass to fn as arguments
#'@param ... Additional parameters to be passed to fn
#'@param smooth Whether to return a smoothed distribution over continuous parameters
#'
#'@export
setGeneric("extract.simset.distribution",
           function(simset,
                    fn,
                    pass.to.fn=c('simulation','parameters','both')[1],
                    ...,
                    smooth=F)
           {
               standardGeneric("extract.simset.distribution")
           })
setMethod("extract.simset.distribution",
          signature(simset='simset', fn='function'),
          def=function(simset,
                       fn,
                       pass.to.fn=c('simulation','parameters','both')[1],
                       ...){

              #check pass.to.fn

              if (pass.to.fn == 'simulation')
                  points = sapply(simset@simulations, fn, ...)
              else if (pass.to.fn == 'parameters')
                  points = apply(simset@parameters, 1, fn, ...)
              else
                  points = sapply(1:simset@n.sim, function(i){
                      fn(simset@simulations[[i]], simset@parameters[i,], ...)
                  })

              if (is.null(dim(points)))
                  points = matrix(points, ncol=1)
              else
                  points = t(points)

              if (!smooth)
                  distributions::Empiric.Distribution(points, weights=simset@weights)
              else
                  stop("Have not yet implemented smoothed distributions. Coming soon")
          })

#'@title Extract the distribution of parameters in a simset
#'
#'@description Gets a probability distribution over the parameters in a simset
#'
#'@inheritParams extract.simset.distribution
#'@param simset A \link{simset} object
#'@param smooth Whether to return a smoothed distribution over continuous parameters
#'
#'@export
setGeneric("extract.simset.parameter.distribution",
           function(simset, smooth=F)
           {
               standardGeneric("extract.simset.parameter.distribution")
           })
setMethod("extract.simset.parameter.distribution",
          signature(simset='simset'),
          def=function(simset, smooth=F)
          {
              if (!smooth)
                  distributions::Empiric.Distribution(simset@parameters, weights=simset@weights)
              else
                  stop("Have not yet implemented smoothed distributions. Coming soon")
          })

#'@title Flatten a simset
#'
#'@description Flatten a simset such that the weight of each simulation is 1. Simulations in the input simset with weight >1 are repeated multiple times in the returned simset
#'
#'@param simset An object of class simset. Must have ony integer-valued weights
#'
#'@return A new simset object
#'
#'@export
setGeneric("flatten.simset",
           def=function(simset){
               standardGeneric("flatten.simset")
           })
setMethod("flatten.simset",
          signature(simset="simset"),
def=function(simset){
    if (any(round(simset@weights) != simset@weights))
        stop("In order to flatten the simset, all weights must be integers")

    new.indices = unlist(sapply(1:simset@n.sim, function(i){
        rep(i, simset@weights[i])
    }))

    new("simset",
        simulations = simset@simulations[new.indices],
        parameter.names=simset@parameter.names,
        n.parameters=simset@n.parameters,
        parameter.values = simset@parameters[new.indices,],
        n.sim = length(new.indices),
        weights = rep(1, length(new.indices)),
        parameter.transformations = simset@parameter.transformations#,
#        parameter.lower.bounds = simset@parameter.lower.bounds,
 #       parameter.upper.bounds = simset@parameter.upper.bounds
        )
})

##------------------------##
##-- SUBSETTING SIMSETS --##
##------------------------##

#'@title Subset a simset
#'
#'@description Subset a simset by indexing a subset of the included simulations
#'
#'@param simset An object of class simset
#'@param indices Either numeric indices or a logical mask according to which to subset simulations
#'
#'@return A new simset object
#'
#'@export
setGeneric("subset.simset",
           def=function(simset, indices){
               standardGeneric("subset.simset")
           })
setMethod("subset.simset",
          signature(simset="simset"),
def=function(simset, indices){
    simset@simulations = simset@simulations[indices]
    simset@n.sim = length(simset@simulations)
    simset@weights = simset@weights[indices]
    simset@parameters = simset@parameters[indices,]

    simset
})

#'@title Subset a simset by thinning its simulations
#'
#'@description For a simset whose weights are integer values, extracts a subset of the simualtions by thinning every nth simulation (eg, if thin=3, the resulting simset keeps every third simulation). Maintains the weighting - simulations with higher weights are more likely to be kept
#'
#'@param simset A simset object
#'@param thin How much to thin the simulations
#'
#'@return A thinned simset object
#'
#'@export
setGeneric("thin.simset",
           def=function(simset, thin){
               standardGeneric("thin.simset")
           })
setMethod("thin.simset",
          signature(simset="simset"),
def=function(simset, thin)
{
    if (any(round(simset@weights) != simset@weights))
        stop("thin.simset can only be called when all weights are integers")

    if (thin <= 0)
        stop("thin must be >= 1")
    if (thin >sum(simset@weights))
        stop(paste0("A thin of ", thin, " will leave no simulations in the thinned simset"))

    expanded.indices = unlist(sapply(1:simset@n.sim, function(i){
        rep(i, simset@weights[i])
    }))

    keep.expanded.indices.mask = ((1:length(expanded.indices)) %% thin) == 0
    keep.indices = unique(expanded.indices[keep.expanded.indices.mask])
    keep.index.weights = as.numeric(table(expanded.indices[keep.expanded.indices.mask]))

    simset = subset.simset(simset, keep.indices)
    simset@weights = keep.index.weights

    simset
})
##-----------------------##
##-- EXTENDING SIMSETS --##
##-----------------------##

#'@title Add parameters to a simset or parameter_set
#'
#'@param simset A \link{simset} or \link{parameter_set} to which to add parameters
#'@param parameters The parameters to add. If only one new parameter, can be either a vector of length simset@n.sim or length 1, or a matrix with one column and simset@n.sim rows. If more than one parameter is to be added, a matrix with one column for parameter and simset@n.sim rows
#'@param parameter.names The names of the new parameters being added. If NULL, takes the parameter names from the column names of parameters
#'@param parameter.lower.bounds,parameter.upper.bounds The min and max values each parameter can take. If given a scalar value, will apply the same bounds to all parameters
#'@param parameter.transformations Objects of class \link{distributions::transformation}, or the names of defined transformations, according to which parameters should be smoothed when making smoothed distributions
#'
#'@export
setGeneric('add.parameters',
           def=function(simset, parameters, parameter.names=NULL,
                        parameter.transformations=NULL,
                        parameter.lower.bounds=-Inf,
                        parameter.upper.bounds=Inf)
               {standardGeneric('add.parameters')})

setMethod('add.parameters',
          signature(simset='parameter_set'),
def=function(simset, parameters, parameter.names=NULL,
             parameter.transformations=NULL,
             parameter.lower.bounds=-Inf,
             parameter.upper.bounds=Inf)
{
    if (is.null(dim(parameters)))
    {
        if (length(parameters)==1)
            parameters = rep(parameters, simset@n.sim)

        parameters = matrix(parameters, ncol=1)
    }

    if (!is(parameters, 'matrix') && !is(parameters, 'array') &&
        length(dim(parameters)) != 2)
        stop("'parameters' must be either a vector, matrix, or 2-dimensional array")

    if (dim(parameters)[1] != simset@n.sim)
        stop(paste0("There must be simset@n.sim (",
                    simset@n.sim, ") rows/values of parameters"))

    if (is.null(parameter.names))
        parameter.names = dimnames(parameters)[[2]]

    if (!is.null(parameter.names) && length(parameter.names) != dim(parameters)[2])
        stop("If parameter.names is specified, there must be one value for each new parameter (ie one value for each column in parameters)")

    if (is.null(parameter.names) && !is.null(simset@parameter.names))
        stop(paste0("The given simset has named parameters, so additional parameters must have names as well"))

    else if (!is.null(parameter.names) && is.null(simset@parameter.names))
    {
        parameter.names = NULL
        warning("The given simset does not have named parameters; parameter names for the new parameters will be ignored")
    }

    if (is.null(parameter.names))
        dimnames(parameters) = NULL
    else
        dimnames(parameters) = list(NULL, parameter.names)

    simset@parameters = cbind(simset@parameters, parameters)
    simset@parameter.names = c(simset@parameter.names, parameter.names)
    if (!is.null(simset@parameter.names) && max(table(simset@parameter.names))>1)
        stop("parameter.names must be unique and distinct from the previous simset@parameter.names")
    simset@n.parameters = dim(simset@parameters)[2]

    #need to do bounds

    simset@parameter.transformations = c(simset@parameter.transformations,
                                         process.transformations(parameter.names,
                                                                 parameter.transformations,
                                                                 var.names.name='parameter.names'))

    simset
})


#'@title Extend the simulations of a simset
#'
#'@param simset An object of class \link{simset}
#'@param fn A function that takes two arguments: sim (the prior simulation object) and parameters (the corresponding numeric vector of parameters) and returns a new simulation object
#'
#'@return An updated simset
#'
#'@export
setGeneric('extend.simulations',
           def=function(simset, fn){standardGeneric('extend.simulations')})
setMethod('extend.simulations',
          signature(simset='simset'),
def=function(simset, fn){
    simset@simulations = lapply(1:simset@n.sim, function(i){
        fn(simset@simulations[[i]], simset@parameters[i,])
    })

    simset
})
