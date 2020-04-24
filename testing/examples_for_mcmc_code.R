

## To run from a cache ##

#This runs it with a cache. Same arguments as run.mcmc, with three extras
# cache.frequency, cache.dir, remove.cache.when.done

mcmc = run.mcmc.(ctrl,
                 n.iter=N.ITER,
                 starting.values=starting.values,
                 update.frequency=1000,
                 cache.dir = 'wherever you want to temporarily cache files',
                 cache.frequency = 500,
                 remove.cache.when.done = T)

mcmc = assemble.mcmc.from.cache(dir='wherever you want to temporarily cache files', allow.incomplete = T)



##-- To analyze the MCMC chain --##

library(mcmc.sim)
#I took a partially run mcmc from my hiv model for you
# (I simplified the simulations objects so as not to take up a lot of file space)
load('example.mcmc.Rdata')

mcmc@n.iter #see how many iterations we have here (after thinning)

#-- Trace plots - the trace.plot function --#
trace.plot(mcmc) #by default, this plots all the variables. That's a lot here

trace.plot(mcmc, c('peak.hiv.mortality','hiv.mortality.0','hiv.mortality.1','hiv.mortality.2'))
        #or I can name which specific variables to plot - here, the four mortality parameters (a spline)

trace.plot(mcmc, '*mortality')
    #but that always seems like a lot of typing (and to remember the exact names)
    #so you can also just do wild card matching

#if you prefer to access the sampled values directly, you can use
mcmc@samples[chain, iteration, variable.name]


#if you have more than one chain, you can calculate Gelmna's rhat statistic to monitor convergence
get.rhats(mcmc)

#-- Log likelihood and prior --#

#we store the log likelihood and log prior, accessed [chain,iteration]
mcmc@log.likelihoods[1,]
mcmc@log.priors[1,]

#or we can plot them
likelihood.plot(mcmc)

#since the prior is often quite different from the likelihood, it's useful to plot them separately
likelihood.plot(mcmc, show.log.prior = F)
likelihood.plot(mcmc, show.log.likelihood = F, show.log.prior.plus.likelihood = F)

#-- Acceptance rates --#

#to check how our acceptance rate over the past 50 iterations
get.cumulative.acceptance.rates(mcmc, window.iterations=50)

#or you can just plot it
acceptance.plot(mcmc, window.iterations = 50)


#-- Turn it into a simset --#
simset = extract.simset(mcmc)

#if you want to burn extra initial simulations (ie, not use the first simulations), use additional.burn
#If you want to thin further, use additional.thin
simset = extract.simset(mcmc, additional.burn=100, additional.thin=2)

#we can get a posterior distribution over our parameters, and do inference from it
post.param.dist = extract.simset.parameter.distribution(simset)
get.means(post.param.dist)
get.intervals(post.param.dist, coverage=.95)

#we can do the same sort of thing - get posterior distributions - for any quantity
# we extract from our simulation objects (returned by run.simulation)
#for example, I'm going to extract the notifications from 2010-2020 from my simplified simulation objects:
get.sim.notifications <- function(sim)
{
    rv = sim[as.character(2010:2020)]
    names(rv) = paste0('notifications_', 2010:2020)
    rv
}
notifications.distribution = extract.simset.distribution(simset, fn=get.sim.notifications)
notifications.distribution@var.names
get.means(notifications.distribution)
get.intervals(notifications.distribution)
