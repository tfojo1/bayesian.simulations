
#'@title Generate a trace plot from an mcmc run
#'
#'@param mcmc A \link{mcmc.sim} object
#'@param var.names,exact.var.names var.names is the names of variables to include or a set of integer indices. If exact.var.names is TRUE, matches the given strings exactly. If exact.var.names is FALSE, matches any variable name that begins with any element of var.names, where the * character is treated as a wild card
#'@param chains Which chains from the MCMC to use
#'@param ncol Number of columns in the facet_wrap for the resulting plot object
#'@param include.rhat If TRUE, includes the rhat value in the facet header
#'
#'@export
trace.plot <- function(mcmc,
                       var.names=mcmc@var.names,
                       chains = 1:mcmc@n.chains,
                       additional.burn=0,
                       additional.thin=1,
                       exact.var.names=F,
                       include.rhat=mcmc@n.chains>1,
                       ncol=NULL)
{
    #Pull samples
    samples = do.subset.mcmc(mcmc, chains=chains, var.names=var.names,
                                   additional.burn=additional.burn, additional.thin=additional.thin)$samples

    if (include.rhat)
    {
        rhats = get.rhats(mcmc, var.names=dimnames(samples)[['variable']],
                          chains=chains, additional.burn=additional.burn,
                          additional.thin=additional.thin, exact.var.names = T)[dimnames(samples)[['variable']]]

        dimnames(samples)[['variable']] = paste0(dimnames(samples)[['variable']],
                                                 " (",
                                                 format(round(rhats, 2), nsmall = 2),
                                                 ")")
    }
    #Make our data object
#    samples = pared.samples@samples[,,var.names]
#    dim(samples) = c(chain=length(chains), iteration=mcmc@n.iter, variable=length(var.names))
#    dimnames(samples) = list(chain=as.character(chains), iteration=NULL, variable=var.names)

    df = reshape2::melt(samples)
    df$chain = factor(df$chain, levels=chains)
    df$iteration = as.integer(df$iteration)
    ggplot2::ggplot(df, ggplot2::aes(iteration, value, color=chain)) +
        ggplot2::geom_line() + ggplot2::facet_wrap(~variable, scales='free_y', ncol=ncol)
}

#'@title Generate a pairs plot of variables from an mcmc run
#'
#'@inheritParams trace.plot
#'
#'@export
pair.plot <- function(mcmc,
                       var.names=mcmc@var.names,
                       chains = 1:mcmc@n.chains,
                       additional.burn=0,
                       additional.thin=1,
                       exact.var.names=F)
{
    #Pull samples
    samples = do.subset.mcmc(mcmc, chains=chains, var.names=var.names,
                             additional.burn=additional.burn, additional.thin=additional.thin)$samples

    df = NULL
    for (var1 in dimnames(samples)[['variable']])
    {
        for (var2 in dimnames(samples)[['variable']])
        {
            df1 = reshape2::melt(samples[,,var1])
            df2 = reshape2::melt(samples[,,var2])
            df = rbind(df,
                       data.frame(value1 = df1$value,
                                  value2 = df2$value,
                                  chain = df1$chain,
                                  var1=var1,
                                  var2=var2
                                  ))
        }
    }
    n.var = dim(samples)['variable']

    df$chain = factor(df$chain, levels=chains)
    ggplot2::ggplot(df, ggplot2::aes(value1, value2, color=chain)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~var1+var2, scales='free', ncol=n.var)
}

density.plot <- function(mcmc)
{

}

#'@export
acceptance.plot <- function(mcmc,
                            window.iterations=ceiling(200/mcmc@thin),
                            by.block=F,
                            aggregate.chains=F,
                            blocks=mcmc@sample.steps,
                            exact.block.names=F)
{
    if (by.block)
    {
        if (!aggregate.chains)
            warning("Cannot show chains and block simultaneously. Showing blocks only")
        aggregate.chains=T
    }

    rates = get.cumulative.acceptance.rates(mcmc=mcmc,
                                       window.iterations=window.iterations,
                                       by.block=by.block,
                                       aggregate.chains=aggregate.chains,
                                       blocks=blocks,
                                       exact.block.names=exact.block.names)
    df = reshape2::melt(rates)
    df = df[!is.na(df$value),]

    if (!aggregate.chains)
        df$chain = factor(df$chain)

    if (window.iterations >= mcmc@n.iter)
        y.lab = "Cumulative Acceptance Rate"
    else
        y.lab = paste0("Acceptance Rate over prior ",
                       format(window.iterations, big.mark=','),
                       " iterations")

    if (by.block)
        rv = ggplot2::ggplot(df, ggplot2::aes(iteration, value, color=block))
    else if (!aggregate.chains)
        rv = ggplot2::ggplot(df, ggplot2::aes(iteration, value, color=chain))
    else
        rv = ggplot2::ggplot(df, ggplot2::aes(iteration, value))

#    if (.hasSlot(mcmc@control, 'target.acceptance.probability'))
#        rv = rv + ggplot2::geom_hline(yintercept = mcmc@control@target.acceptance.probability, linetype='dashed')

    rv = rv + ggplot2::geom_line()

    if (by.block && !aggregate.chains && mcmc@n.chains>1)
        rv = rv + ggplot2::facet_wrap(~chain)

    rv = rv + ggplot2::ylim(0,1) + ggplot2::ylab(y.lab)

    rv
}

#'@export
likelihood.plot <- function(mcmc,
                            chains = 1:mcmc@n.chains,
                            show.log.prior=T,
                            show.log.likelihood=T,
                            show.log.prior.plus.likelihood=T)
{
    if (!show.log.likelihood && !show.log.likelihood && !show.log.prior)
        stop("You must set at least one of show.log.prior, show.log.likelihood, or show.log.prior.plus.likelihood to TRUE")

    #Check chains
    chains = check.chains(mcmc, chains)

    df = NULL
    if (show.log.likelihood)
    {
        one.df = reshape2::melt(mcmc@log.likelihoods)
        one.df$type='Log Likelihood'

        df = rbind(df, one.df)
    }

    if (show.log.prior)
    {
        one.df = reshape2::melt(mcmc@log.priors)
        one.df$type='Log Prior'

        df = rbind(df, one.df)
    }

    if (show.log.prior.plus.likelihood)
    {
        one.df = reshape2::melt(mcmc@log.likelihoods + mcmc@log.priors)
        one.df$type='Log Prior x Likelihood'

        df = rbind(df, one.df)
    }

    #select chains
    df = df[sapply(df$chain, function(ch){any(ch==chains)}),]
    df$chain = factor(df$chain, levels=sort(unique(df$chain)))

    ggplot2::ggplot(df, ggplot2::aes(iteration, value, color=chain, linetype=type)) +
        ggplot2::geom_line()
}

#'@export
#'
#'@return Returns an array indexed [chain,iteration,block]. If by.block==T, then there is one value of block for each var.block. If by.block==F, then there is only a single value of block ('all'). The values of the array represent the cumulative fraction of transitions that have been accepted from the end of the burn period up to and including the current iteration (including thinned transitions)
get.cumulative.acceptance.rates <- function(mcmc,
                                               window.iterations=ceiling(200/mcmc@thin),
                                               by.block=F,
                                               aggregate.chains=F,
                                               blocks=mcmc@sample.steps,
                                               exact.block.names=F)
{
    blocks = match.block.names(mcmc, blocks, exact.block.names)
    dim.names = list(chain=1:mcmc@n.chains, iteration=1:mcmc@n.iter, block=blocks)

    n.accepted = mcmc@n.accepted[,,blocks]
    dim(n.accepted) = sapply(dim.names, length)
    dimnames(n.accepted) = dim.names

    cum.n.accepted = apply(n.accepted, c(1,3), cumsum)
    cum.n.accepted = apply(cum.n.accepted, c(2,1,3), function(x){x})
    dim(cum.n.accepted) = sapply(dim.names, length)
    dimnames(cum.n.accepted) = dim.names

    n.unthinned = get.n.unthinned.by.block(mcmc)[,,blocks]
    dim(n.unthinned) = sapply(dim.names, length)
    dimnames(n.unthinned) = dim.names

    cum.n.unthinned = apply(n.unthinned, c(1,3), cumsum)
    cum.n.unthinned = apply(cum.n.unthinned, c(2,1,3), function(x){x})
    dim(cum.n.unthinned) = sapply(dim.names, length)
    dimnames(cum.n.unthinned) = dim.names

    if (window.iterations < mcmc@n.iter)
    {
        num.to.subtract = mcmc@n.iter-window.iterations
        cum.n.accepted[,window.iterations + 1:num.to.subtract,] = cum.n.accepted[,window.iterations + 1:num.to.subtract,] -
            cum.n.accepted[,1:num.to.subtract,]
        cum.n.unthinned[,window.iterations + 1:num.to.subtract,] = cum.n.unthinned[,window.iterations + 1:num.to.subtract,] -
            cum.n.unthinned[,1:num.to.subtract,]
    }

    rv = apply(cum.n.accepted, (1:3)[c(!aggregate.chains, T, by.block)], sum) /
            apply(cum.n.unthinned, (1:3)[c(!aggregate.chains, T, by.block)], sum)

    dim.names = dim.names[c(!aggregate.chains, T, by.block)]
    dim(rv) = sapply(dim.names, length)
    dimnames(rv) = dim.names

    rv
}

#'@export
get.total.acceptance.rate <- function(mcmc,
                                  by.block=F,
                                  aggregate.chains=T,
                                  blocks=mcmc@sample.steps,
                                  exact.block.names=F)
{
    blocks = match.block.names(mcmc, blocks, exact.block.names)
    dim.names = list(chain=1:mcmc@n.chains, iteration=1:mcmc@n.iter, block=blocks)

    n.accepted = mcmc@n.accepted[,,blocks]
    dim(n.accepted) = sapply(dim.names, length)
    dimnames(n.accepted) = dim.names

    n.unthinned = get.n.unthinned.by.block(mcmc)[,,blocks]
    dim(n.unthinned) = sapply(dim.names, length)
    dimnames(n.unthinned) = dim.names

    if (aggregate.chains && !by.block)
        rv = sum(n.accepted) / sum(n.unthinned)
    else
    {
        rv = apply(n.accepted, (1:3)[c(!aggregate.chains, F, by.block)], sum) /
            apply(n.unthinned, (1:3)[c(!aggregate.chains, F, by.block)], sum)

        dim.names = dim.names[c(!aggregate.chains, F, by.block)]
        dim(rv) = sapply(dim.names, length)
        dimnames(rv) = dim.names
    }

    rv
}

#'@title Calculate R-hat statistics for an MCMC object
#'
#'@description Calculates the R-hat of Gelman 2013, or the R-hat using rank-normalization as per Vehtari 2019
#'
#'@inheritParams trace.plot
#'@param rank.normalize Whether to rank-normalize samples before calculating the R-hat, as in Vehtari 2019
#'@param sort If TRUE, the rhats are sorted from greatest to least
#'
#'@export
get.rhats <- function(mcmc,
                      var.names=mcmc@var.names,
                      chains=1:mcmc@n.chains,
                      additional.burn=0,
                      additional.thin=1,
                      exact.var.names=F,
                      rank.normalize=F,
                      sort=T)
{
    #Check chains
    chains = check.chains(mcmc, chains)
    if (length(chains)==1)
        stop('R-hats can only be calculated with more than one chain')



    #Pull samples
    samples = do.subset.mcmc(mcmc, chains=chains, var.names=var.names,
                             additional.burn=additional.burn, additional.thin=additional.thin)$samples


    samples = apply(samples, c('variable','chain','iteration'), function(x){x})
        #now indexed [variable, chain, iteration]
    dim.names = dimnames(samples)
    var.names = dim.names[['variable']]

    #?rank normalize
    if (rank.normalize)
    {
        samples = apply(samples, 'variable', rank.normalize)
        dim(samples) = sapply(dim.names, length)
        dimnames(samples) = dim.names
    }

    #calculate rhats
    N = mcmc@n.iter
    M = length(chains)

    theta.hat.dot.m = rowMeans(samples, dims=2)
    dim(theta.hat.dot.m) = sapply(dim.names[c('variable','chain')], length)
    dimnames(theta.hat.dot.m) = dim.names[c('variable','chain')]
        #indexed [variable, chain]

    theta.hat.dot.dot = rowMeans(theta.hat.dot.m)
        # indexed [variable]

    B = N / (M-1) * rowSums( (theta.hat.dot.m - theta.hat.dot.dot)^2 )
        # indexed [variable]

    #we're going to flatten samples and theta.hat.dot.m so we can subtract them
    theta.hat.dot.m = as.numeric(theta.hat.dot.m)
    dim(samples) = c(prod(dim(samples)[1:2]), dim(samples)[3])

    s.m.squared = 1 / (N-1) * rowSums( (samples - theta.hat.dot.m)^2 )
    dim(s.m.squared) = sapply(dim.names[c('variable','chain')], length)
    dimnames(s.m.squared) = dim.names[c('variable','chain')]
        # indexed [variable, chain]

    W = rowMeans(s.m.squared)
        # indexed [variable]

    var.hat.plus.theta.given.y = (N-1)/N * W + 1/N * B

    rhats = sqrt(var.hat.plus.theta.given.y/ W)
    names(rhats) = var.names

    if (sort)
        rhats = sort(rhats, decreasing=T)

    rhats
}

#Think about Seff
#https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html

##-- HELPERS --##
get.n.unthinned.by.block <- function(mcmc)
{
    dim.names = list(chain=1:mcmc@n.chains, iteration=1:mcmc@n.iter, block=mcmc@sample.steps)
    rv = array(0, dim=sapply(dim.names, length), dimnames=dim.names)
    n.blocks = length(mcmc@sample.steps)

    for (chain in 1:mcmc@n.chains)
    {
        counts = t(sapply(1:mcmc@n.iter, function(iter){
            raw.counts = 1 + ((mcmc@first.step.for.iter[iter] - 2 + 1:mcmc@thin) %% n.blocks)
            sapply(1:n.blocks, function(block){
                sum(raw.counts==block)
            })
        }))

        rv[chain,,] = counts
    }

    rv
}

check.chains <- function(mcmc, chains)
{
    if (!is(chains, 'integer') &&
        !(is(chains, 'numeric') && all(round(chains)==chains)))
        stop(paste0("'chains' must be set of non-NA integers between 1 and ", mcmc@n.chains))

    if (any(is.na(chains) || any(chains<1) || any(chains>mcmc@n.chains)))
        stop(paste0("'chains' must be set of non-NA integers between 1 and ", mcmc@n.chains))

    unique(chains)
}

rank.normalize <- function(x)
{
    r = rank(x, ties.method = 'average')
    z = qnorm( (r - 3/8) / (length(s) - 1/4))

    if (!is.null(dim(x)))
    {
        dim(z) = dim(x)
        dimnames(z) = dimnames(x)
    }

    z
}
