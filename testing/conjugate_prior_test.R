
library(mcmc.sim)
library(distributions)
library(mvtnorm)

set.seed(1234)

N.ITER = 100000
N.BURN = 50000

#-- SET UP THE TRUTH --#
#   An MVN with mean vector true.means and covariance matrix sigma
n.dim = 40
var.names = paste0('mu', 1:n.dim)
var.blocks = list(c('mu1','mu2'), 'mu3', var.names[4:n.dim])
var.blocks = list(paste0('mu',1:5),
                  paste0('mu',6:7),
                  paste0('mu',8:10),
                  paste0('mu',11:15),
                  paste0('mu',16:17),
                  paste0('mu',18:20),
                  paste0('mu',21:25),
                  paste0('mu',26:27),
                  paste0('mu',28:30),
                  paste0('mu',31:35),
                  paste0('mu',36:37),
                  paste0('mu',38:40))
#var.blocks = as.list(var.names)

N.ITER.BEFORE.COV = 0

rho=0.8
ar.cor = rho^matrix(abs(rep(1:n.dim, each=n.dim)-rep(1:n.dim, n.dim)), nrow=n.dim, ncol=n.dim)

true.means = rnorm(n.dim, 0, 40)
names(true.means) = var.names
sds = runif(n.dim, 1, 20)
sigma = sds %*% t(sds) * ar.cor

#-- SET UP THE DATA --#
n.data = 100
data = rmvnorm(n.data, mean=true.means, sigma = sigma)

#-- SET UP THE PRIOR AND LIKELIHOOD --#
mu0 = rep(0, n.dim)
sigma0 = diag(rep(100, n.dim))
log.prior <- function(mu){
    dmvnorm(mu, mean=mu0, sigma=sigma0, log=T)
}

log.likelihood <- function(mu){
    sum(dmvnorm(data, mean=mu, sigma=sigma, log=T))
}

simulation.function <- function(mu){mu}

#-- SET UP THE STARTING VALUES AND MCMC CONTROL --#
starting.values = rep(0, n.dim)
names(starting.values) = var.names

ctrl = create.adaptive.blockwise.metropolis.control(var.names=var.names,
                                                    var.blocks=var.blocks,
                                                    log.likelihood = log.likelihood,
                                                    log.prior.distribution = log.prior,
                                                    simulation.function = simulation.function,
                                                    initial.covariance.mat = diag(rep(1,n.dim)),
                                                    transformations = NULL,
                                                    n.iter.before.use.adaptive.covariance = N.ITER.BEFORE.COV,
                                                    thin=50,
                                                    burn=N.BURN)
#-- RUN THE MCMC --#
mcmc = run.mcmc(ctrl, n.iter=N.ITER, starting.values = starting.values, update.detail = 'high')
simset = extract.simset(mcmc)
dist = extract.simset.parameter.distribution(simset)

get.means(dist)
get.covariance.matrix(dist)

#-- THE CLOSED FORM POSTERIOR --#
sample.mean = colMeans(data)
sigma.inverse = solve(sigma)
sigma0.inverse = solve(sigma0)

posterior.cov.mat = solve(sigma0.inverse + n.data * sigma.inverse)
posterior.mean = posterior.cov.mat %*%
    (sigma0.inverse %*% mu0 + n.data * sigma.inverse %*% sample.mean)


#-- COMPARISON --#
cbind(mcmc=get.means(dist), closed_form=posterior.mean)
cbind(mcmc=as.numeric(get.covariance.matrix(dist)),
      closed_form=as.numeric(posterior.cov.mat))

library(ggplot2)
qplot(posterior.mean, get.means(dist)) + ggtitle("Means (true posterior vs mcmc)") + geom_abline(intercept = 0, slope = 1)
qplot(as.numeric(get.covariance.matrix(dist)), as.numeric(posterior.cov.mat)) + ggtitle("Covariance Matrix elements (true posterior vs mcmc)") + geom_abline(intercept = 0, slope = 1)

if (1==2)
{
    qplot(posterior.mean, colMeans(mcmc.one@samples[1,,])) + ggtitle("Means (true posterior vs mcmc)") + geom_abline(intercept = 0, slope = 1)
    qplot(as.numeric(posterior.cov.mat), as.numeric(cov(mcmc.one@samples[1,,]))) + ggtitle("Means (true posterior vs mcmc)") + geom_abline(intercept = 0, slope = 1)
}

#-- ANALYSIS --#
mcmc@chain.states[[1]]@log.scaling.parameters
trace.plot(mcmc)

if (1==2)
{
    cbind(mcmc.ncov.single@chain.states[[1]]@log.scaling.parameters,
          mcmc.single.with.cov@chain.states[[1]]@log.scaling.parameters,
          mcmc.with.cov.multi.step@chain.states[[1]]@log.scaling.parameters)
}
