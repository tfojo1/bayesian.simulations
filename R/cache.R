
##------------------------##
##-- THE MASTER CLASSES --##
##------------------------##

setClassUnion('mcmcsim_or_NULL', members = c('NULL', 'mcmcsim'))

setClass('mcmcsim_cache_global_control',
         representation=list(
             id='character',
             control='mcmcsim_control',
             n.chains='integer',
             chain.control.filenames='character',
             n.chunks='integer',
             chunk.size='integer',
             save.chunk='logical',
             need.to.remove.dir='logical',
             prior.mcmc='mcmcsim_or_NULL'
         ))

setClass('mcmcsim_cache_chain_control',
         representation=list(
             global.id='character',
             chain.id='character',
             chunk.filenames='character',
             control.filename='character',
             chain.dir='character',
             chunk.done='logical',
             n.accepted.so.far='integer',
             total.runtime='numeric',
             chain.state='mcmcsim_chainstate'
         ))

##--------------------##
##-- SOME CONSTANTS --##
##--------------------##

CACHE.SUBDIR.NAME = 'cache'
CACHE.GLOBAL.CONTROL.NAME = 'global_control.Rdata'

get.cache.subdirectory <- function(dir)
{
    file.path(dir, CACHE.SUBDIR.NAME)
}

get.cache.global.control.filename <- function(dir)
{
    file.path(dir, CACHE.SUBDIR.NAME, CACHE.GLOBAL.CONTROL.NAME)
}

get.cache.chain.control.filename <- function(dir, global.control, chain)
{
    file.path(get.cache.subdirectory(dir), global.control@chain.control.filenames[chain])
}

##---------------------------------##
##-- CONSISTENCY CHECK FUNCTIONS --##
##---------------------------------##

cache.exists <- function(dir)
{
    dir.exists(get.cache.subdirectory(dir))
#    filename = get.cache.global.control.filename(dir)
#    dir.exists(dir) && file.exists(filename) && is(get(load(filename)), 'mcmcsim_cache_global_control')
}

check.cache.corrupted <- function(dir, chains=NULL)
{
    global.control = get.cache.global.control(dir)
    if (is.null(chains))
        chains = 1:global.control@n.chains
    else if (length(setdiff(chains, 1:global.control@n.chains))>0)
        stop(paste0("chains must be between 1 and ", global.control@n.chains))

    for (chain in chains)
    {
        chain.control = get.cache.chain.control(dir, global.control, chain)

        if (chain.control@global.id != global.control@id)
            stop(paste0("The cache has been corrupted. The control for the ",
                        get.ordinal(chain),
                        " chain does not derive from the global control"))

        chain.dir = file.path(get.cache.subdirectory(dir), chain.control@chain.dir)
        if (!dir.exists(chain.dir))
            stop(paste0("The cache has been corrupted. The directory for the ",
                        get.ordinal(chain),
                        " ('",
                        chain.dir,
                        "') does not exist"))

        chunk.files = file.path(get.cache.subdirectory(dir), chain.control@chain.dir, chain.control@chunk.filenames)
        missing.chunks = (1:global.control@n.chunks)[chain.control@chunk.done &
                                                         global.control@save.chunk &
                                                         !file.exists(chunk.files)]
        if (length(missing.chunks)>0)
            stop(paste0("The cache has been corrupted. The ",
                        paste0(get.ordinal(missing.chunks), collapse=', '),
                        " sampling chunks from chain ",
                        chain,
                        " have been done, but the files containing the cached results are missing: ",
                        paste0("'", chunk.files[missing.chunks], "'", collapse=', ')))

#        extra.chunks = (1:global.control@n.chunks)[!chain.control@chunk.done &
#                                                       global.control@save.chunk &
#                                                       file.exists(chunk.files)]

#        if (length(extra.chunks)>0)
#            stop(paste0("The cache has been corrupted. The ",
#                        paste0(get.ordinal(extra.chunks), collapse=', '),
#                        " sampling chunks from chain ",
#                        chain,
#                        " have NOT been done, but the files to contain their results already exist: ",
#                        paste0("'", chunk.files[extra.chunks], "'", collapse=', ')))

    }
}

##-------------------##
##-- BASIC GETTERS --##
##-------------------##

get.cache.global.control <- function(dir, throw.error.if.not.found=T)
{
    if (!dir.exists(dir))
    {
        if(throw.error.if.not.found)
            stop(paste0("The directory '", dir, '" does not exist'))
        else
            return (NULL)
    }

    filename = get.cache.global.control.filename(dir)
    if (!file.exists(filename))
    {
        if (throw.error.if.not.found)
            stop(paste0("No cache is set up in directory '", dir, "'"))
        else
            return (NULL)
    }

    global.control=get(load(filename))

    if (!is(global.control, 'mcmcsim_cache_global_control'))
    {
        if (throw.error.if.not.found)
            stop(paste0("'", filename, "' does not contain an mcmcsim_cache_global_control object"))
        else
            return (NULL)
    }

    global.control
}

get.cache.chain.control <- function(dir, global.control, chain)
{
    if (!is(global.control, 'mcmcsim_cache_global_control'))
        stop("global.control must be an object of class 'mcmcsim_cache_global_control'")

    if ((!is(chain, 'integer') && !is(chain, 'numeric')) || length(chain)!=1)
        stop("chain must be a scalar integer value")

    if (chain < 1 || chain > global.control@n.chains)
        stop(paste0('chain must be between 1 and ', global.control@n.chains))

    filename = get.cache.chain.control.filename(dir, global.control, chain)

    if (!file.exists(filename))
        stop(paste0("The control file for the ", get.ordinal(chain), " chain is missing"))

    chain.control = get(load(filename))

    if (!is(chain.control, 'mcmcsim_cache_chain_control'))
        stop(paste0("'", filename, "' does not contain an mcmcsim_cache_chain_control object"))

    chain.control
}

##--------------------##
##-- CREATE A CACHE --##
##--------------------##

do.create.cache <- function(dir,
                            control,
                            chain.states,
                            n.iter,
                            cache.frequency,
                            prior.mcmc)
{
    #-- Check Cache Frequency --#
    if (is.null(cache.frequency) ||
        (!is(cache.frequency, 'numeric') && !is(cache.frequency('integer'))) ||
        is.na(cache.frequency) ||
        length(cache.frequency) != 1 ||
        cache.frequency%%1 != 0)
        stop("cache.frequency must be a scalar integer")

    #-- Set up chunks and sizes --#
    n.chunks = ceiling(n.iter/cache.frequency)
    chunk.size = c(rep(cache.frequency, n.chunks-1), n.iter-cache.frequency*(n.chunks-1))

    need.to.save.chunk = cumsum(chunk.size) > control@burn
    num.to.save = sum(need.to.save.chunk)

    #-- Create Global Control Object --#
    n.chains = length(chain.states)
    global.control = new('mcmcsim_cache_global_control',
                         id=paste0('global_', Sys.time()),
                         control=control,
                         n.chains=length(chain.states),
                         chain.control.filenames=paste0('chain', 1:n.chains, '_control.Rdata'),
                         chunk.size=as.integer(chunk.size),
                         n.chunks=as.integer(n.chunks),
                         save.chunk=need.to.save.chunk,
                         need.to.remove.dir=F,
                         prior.mcmc=prior.mcmc
                         )

    #-- Set up the Chain Control Objects --#
    chain.controls = do.create.chain.cache.controls(global.control, chains=1:global.control@n.chains,
                                   chain.states=chain.states)

    #-- Create the files and directories --#
    do.create.cache.from.controls(dir,
                                  global.control,
                                  chains=1:n.chains,
                                  chain.controls=chain.controls,
                                  make.global=T)
}


do.create.chain.cache.controls <- function(global.control, chains,
                                           chain.states)
{
    num.to.save = sum(global.control@save.chunk)
    lapply(chains, function(chain){

        chain.dir = paste0('chain_', chain)

        chunk.filenames = paste0('chain', chain, "_chunk", 1:global.control@n.chunks, '.Rdata')
        chunk.filenames[!global.control@save.chunk] = NA

        new('mcmcsim_cache_chain_control',
            global.id = global.control@id,
            chain.id = paste0("chain", chain, "_", Sys.time()),
            chunk.filenames=chunk.filenames,
            control.filename=global.control@chain.control.filenames[chain],
            chain.dir=chain.dir,
            chunk.done=rep(F, global.control@n.chunks),
            n.accepted.so.far=as.integer(0),
            total.runtime=0,
            chain.state=chain.states[[chain]])
    })
}

do.create.cache.from.controls <- function(dir,
                                          global.control,
                                          chains,
                                          chain.controls,
                                          make.global)
{
    #-- Set Up the Global Directory --#

    if (make.global)
    {
        if (!dir.exists(dir))
        {
            dir.created = dir.create(dir)
            if (!dir.created)
                stop("The given directory, '", dir, "' does not exist and could not be created.")
        }
        else
        {
            old.global.control = get.cache.global.control(dir, throw.error.if.not.found = F)
            dir.created = !is.null(old.global.control) && old.global.control@need.to.remove.dir
        }

        if (dir.exists(get.cache.subdirectory(dir)))
            stop(paste0("A cache already exists in directory '", dir, "'"))

        if (!dir.create(get.cache.subdirectory(dir)))
            stop(paste0("Unable to create the cache subdirectory at '", get.cache.subdirectory(dir)))

        global.control@need.to.remove.dir = dir.created
        save(global.control, file=get.cache.global.control.filename(dir))
    }

    for (i in 1:length(chain.controls))
    {
        chain.control = chain.controls[[i]]
        chain = chains[i]
        chain.dir = file.path(get.cache.subdirectory(dir), chain.control@chain.dir)
        if (dir.exists(chain.dir))
            stop(paste0("The directory for storing chunks for the ",
                        get.ordinal(chain), " chain ('",
                        chain.dir,
                        "') already exists"))
        if (!dir.create(chain.dir))
            stop(paste0("Unable to create the directory for storing chunks for the ",
                        get.ordinal(chain), " chain ('",
                        chain.dir,
                        "')"))

        save(chain.control, file=file.path(get.cache.subdirectory(dir), chain.control@control.filename))
    }
}

##--------------------##
##-- REMOVE A CACHE --##
##--------------------##

do.remove.cache <- function(dir)
{
    if (dir.exists(get.cache.subdirectory(dir)))
    {
        global.control = get.cache.global.control(dir, throw.error.if.not.found = F)
        if (unlink(get.cache.subdirectory(dir), recursive = T))
            print(paste0("Unable to delete the cache subdirectory at '",
                        get.cache.subdirectory(dir), "'"))
        else if (!is.null(global.control) && global.control@need.to.remove.dir)
        {
            if(length(list.files(dir, all.files=TRUE, no..=T)) == 0)
            {
                if (unlink(dir, recursive=T))
                    print(paste0("Unable to delete directory '",
                                 dir,
                                 "'. It's contents were deleted"))
            }
            else
                print(paste0("Unable to delete directory '",
                             dir,
                             "' because it is not empty even after deleting the cache"))
        }
    }
}

##----------------------------##
##-- RUN A CHAIN WITH CACHE --##
##----------------------------##

do.run.single.chain.with.cache <- function(dir,
                                           global.control,
                                           chain,
                                           update.frequency,
                                           update.detail,
                                           output.stream)
{
    chain.control = get.cache.chain.control(dir, global.control, chain)

    #-- Run the loop --#
    initial.sim = NULL
    total.n.iter = sum(global.control@chunk.size)
    total.n.iter.undone = sum(global.control@chunk.size[!chain.control@chunk.done])
    n.iter.so.far = sum(global.control@chunk.size[chain.control@chunk.done])

    if (!is.na(update.frequency))
        output.stream(paste0("PREPARING TO RUN ",
                     format(total.n.iter.undone, big.mark=','),
                     " ITERATIONS ACROSS ",
                     sum(!chain.control@chunk.done),
                     " CHUNKS\n"))

    start.time = Sys.time()
    while (any(!chain.control@chunk.done))
    {
        chunk = (1:global.control@n.chunks)[!chain.control@chunk.done][1]

        # Run MCMC
        mcmc.and.sim = run.single.chain(control=global.control@control,
                                        chain.state=chain.control@chain.state,
                                        n.iter=global.control@chunk.size[chunk],
                                        update.frequency=update.frequency,
                                        update.detail=update.detail,
                                        initial.sim=initial.sim,
                                        total.n.iter=total.n.iter,
                                        prior.n.iter=n.iter.so.far,
                                        prior.n.accepted=chain.control@n.accepted.so.far,
                                        prior.run.time=chain.control@total.runtime,
                                        return.current.sim=T,
                                        chain=chain,
                                        output.stream=output.stream)

        # Save MCMC to disk
        if (global.control@save.chunk[chunk])
        {
            mcmc = mcmc.and.sim$mcmc
            save(mcmc, file=file.path(get.cache.subdirectory(dir), chain.control@chain.dir, chain.control@chunk.filenames[chunk]))
        }

        # Update the cache control and save to disk
        chain.control@chain.state = mcmc.and.sim$mcmc@chain.states[[1]]
        chain.control@chunk.done[chunk] = T
        chain.control@n.accepted.so.far = chain.control@n.accepted.so.far +
            as.integer(sum(mcmc.and.sim$mcmc@n.accepted) + sum(mcmc.and.sim$mcmc@n.accepted.in.burn))
        chain.control@total.runtime = chain.control@total.runtime +
            as.numeric(difftime(Sys.time(), start.time, units='secs'))
        start.time = Sys.time()

        save(chain.control, file=get.cache.chain.control.filename(dir, global.control, chain))

        # Update other state variables
        initial.sim = mcmc.and.sim$current.sim
        n.iter.so.far = n.iter.so.far + global.control@chunk.size[chunk]

        # Clear the memory devoted to the last MCMC
        mcmc = NULL
        mcmc.and.sim = NULL
        gc()
    }
}

##---------------------------##
##-- ASSEMBLING FROM CACHE --##
##---------------------------##

do.assemble.single.chain.from.cache <- function(dir,
                                                global.control,
                                                chain,
                                                chunks)
{
    chain.control = get.cache.chain.control(dir, global.control, chain)

    chunks = sort(intersect(chunks, (1:global.control@n.chunks)[global.control@save.chunk]))
    mcmc.objects = lapply(chain.control@chunk.filenames[chunks], function(filename){
        load(file.path(get.cache.subdirectory(dir), chain.control@chain.dir, filename))
        mcmc
    })
    mcmc = mcmc.merge.serial(mcmc.objects)


    #-- Return --#
    mcmc
}

do.assemble.multiple.chains.from.cache <- function(dir,
                                                   global.control,
                                                   chains,
                                                   chunks)
{
    single.chain.mcmcs = lapply(chains, function(chain){
        do.assemble.single.chain.from.cache(dir,
                                            global.control=global.control,
                                            chain=chain,
                                            chunks=chunks)
    })

    if (length(single.chain.mcmcs)==1)
        rv = single.chain.mcmcs[[1]]
    else
        rv = mcmc.merge.parallel(single.chain.mcmcs)

    # Merge serial
    if (!is.null(global.control@prior.mcmc))
    {
        #        rv = mcmc.merge.serial(control@prior.mcmc, rv)
        print('need to subset prior mcmc if available')
    }


    rv
}

##---------------------------##
##-- CHECKING CACHE STATUS --##
##---------------------------##

#'@title Check if a cached MCMC has completed sampling
#'
#'@param dir The cache directory
#'
#'@return A boolean indicator of whether the cache is done sampling
#'
#'@export
is.mcmc.cache.complete <- function(dir)
{
    is.cache.complete(dir)
}


is.cache.complete <- function(dir, global.control=NULL, chains=NULL)
{
    # Load cache control
    if (!cache.exists(dir))
        stop("No cache for MCMC exists in directory '", dir, "'")

    if (is.null(global.control))
        global.control = get.cache.global.control(dir)


    # Check chains
    if (is.null(chains))
        chains = 1:global.control@n.chains
    else if (length(setdiff(chains, 1:global.control@chains))>0)
        stop(paste0("chains must include only numbers from 1 to ", global.control@n.chains))

    # Check each chain
    for (chain in chains)
    {
        chain.control = get.cache.chain.control(dir, global.control, chain)
        if (!all(chain.control@chunk.done))
            return (F)
    }

    return (T)
}


get.cache.sampling.status <- function(dir, global.control=NULL, chains=NULL)
{
    # Load cache control
    if (!cache.exists(dir))
        stop("No cache for MCMC exists in directory '", dir, "'")

    if (is.null(global.control))
        global.control = get.cache.global.control(dir)


    # Check chains
    if (is.null(chains))
        chains = 1:global.control@n.chains
    else if (length(setdiff(chains, 1:global.control@chains))>0)
        stop(paste0("chains must include only numbers from 1 to ", global.control@n.chains))

    # Check each chain
    rv = list(chunks.done=rep(0, length(chains)),
              total.chunks=rep(global.control@n.chunks, length(chains)),
              iterations.done=rep(0, length(chains)),
              total.iterations=rep(sum(global.control@chunk.size), length(chains))
              )

    for (i in 1:length(chains))
    {
        chain = chains[i]
        chain.control = get.cache.chain.control(dir, global.control, chain)
        rv$chunks.done[i] = sum(chain.control@chunk.done)
        rv$iterations.done[i] = sum(chain.control@chunk.done * global.control@chunk.size)
    }

    rv
}


chunks.done.per.chain <- function(dir, global.control, chains)
{
    sapply(chains, function(chain){
        chain.control = get.cache.chain.control(dir, global.control, chain)
        cum.chunks.done = cumsum(chain.control@chunk.done)
        chunks.done = (1:global.control@n.chunks)[cum.chunks.done == (1:global.control@n.chunks)]
        max(0, chunks.done)
    })
}
