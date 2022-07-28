
library(ggplot2)

#The two packages I have written
library(bayesian.simulations) #runs the MCMC
library(distributions) # defines distributions as objects - handy for processing results (and potentially defining priors)
                        # not strictly required to run the mcmc

HAVE.NOT.INSTALLED.PACKAGES = F #This is what you need to do to install my packages
if (HAVE.NOT.INSTALLED.PACKAGES)
{
    library(devtools)
    install_github(repo='tfojo1/distributions')
    install_github(repo='tfojo1/bayesian.simulations')
}

##--------------------##
##-- INITIAL SET UP --##
##--------------------##


set.seed = 898798798


##------------------------------------##
##-- SIMULATE 'TRUE' DATA TO FIT TO --##
##------------------------------------##

# Make up some 'true' data we are fitting to
YEARS = 2010:2020
NOTIFICATIONS = rnorm(n=length(YEARS), mean=200 - 5*(YEARS-2010), sd=10)
print(qplot(YEARS, NOTIFICATIONS, geom='line') + ggtitle("Our 'True' Data to fit to"))


##------------------------------------##
##-- DEFINE THE SIMULATION FUNCTION --##
##------------------------------------##

# Defined a toy simulation function
# - take a vector of parameters as input
# - returns a simulation object - however we want to define it
#
# This just projects a line based off two parameters - slope and intercept - and says that the simulated
#  notifications are whatever is produced from that linear model
# In this case, the simulation object is a list where each value is the 'simulated' notifications
#  for a year
run.simulation = function(parameters)
{
    list(simulated.notifications=parameters['intercept'] + parameters['slope'] * (YEARS-2010))
}

##----------------------------##
##-- THE PRIOR DISTRIBUTION --##
##----------------------------##

# Define the log prior for the two parameters
#  intercept ~ normal(250, 50)
#  slope ~ normal(0, 25)
log.prior = function(parameters)
{
    dnorm(parameters['intercept'], mean=250, sd=50, log=T) +
        dnorm(parameters['slope'], mean=0, sd=25, log=T)
}

# This is to show off the distributions package
#  It does the same likelihood but with my nice distribution objects
prior.distribution = join.distributions(
    intercept = Normal.Distribution(mean=250, sd=50),
    slope = Normal.Distribution(mean=0, sd=25)
)

cat("Calculating the likelihood with two versions: \n",
    " Simple prior = ", log.prior(c(intercept=200,slope=-5)),
    " \n  Fancy object prior = ", calculate.density(prior.distribution, c(intercept=200,slope=-5), log=T),
    "\n")

##--------------------##
##-- THE LIKELIHOOD --##
##--------------------##

# Define a simple likelihood
#   For each data point in NOTIFICATIONS, assume the likelihood is
#   NOTIFICATIONS[i] ~ Normal(simulated_notification[i], 25)
log.likelihood = function(sim)
{
    sum(dnorm(NOTIFICATIONS, mean=sim$simulated.notifications, sd=25, log=T))
}


##------------------##
##-- RUN THE MCMC --##
##------------------##


ctrl = create.adaptive.blockwise.metropolis.control(var.names=c('intercept','slope'),
                                                    simulation.function = run.simulation,
                                                    log.prior.distribution = log.prior,
                                                    #To use the distributions package
                                                    #log.prior.distribution = get.density.function(prior.distribution, default.log=T),
                                                    log.likelihood = log.likelihood,
                                                    burn = 1000,
                                                    thin=25,
                                                    initial.covariance.mat = diag(c(10,1)))

#Going to run four chains
n.chains = 4
starting.values = cbind(intercept=rnorm(n.chains, 250, 50),
                        slope=rnorm(n.chains, 0, 25))

RUN.ALL.FOUR.AT.ONCE = F #just to demonstrate you can run these on parallel cores with one call
                         # or separately and then merge later
N.ITER = 5000

CACHE.DIR = '../demonstration_cache'
if (RUN.ALL.FOUR.AT.ONCE) #on parallel cores
{
    mcmc = run.mcmc(ctrl,
                    n.iter=N.ITER,
                    starting.values=starting.values,
                    update.frequency=1000,
                    cache.dir = CACHE.DIR, #this caches as you go, so you can recover if you crash, or 'peek' in
                    cache.frequency=500) #how often to write to a cache
}

if (!RUN.ALL.FOUR.AT.ONCE) #we'll do this using a disk cache
{
    #this saves to a cache
    create.mcmc.cache(dir=CACHE.DIR,
                      control=ctrl,
                      n.iter=N.ITER,
                      starting.values = starting.values,
                      cache.frequency = 500)

    #This is how you could parallelize on a cluster or computer:
    # instead of making each of these calls in serial, you could make one call from each process
    run.mcmc.from.cache(dir=CACHE.DIR,
                        chains=1,
                        update.frequency = 1000) #how often to print updates
    run.mcmc.from.cache(dir=CACHE.DIR,
                        chains=2,
                        update.frequency = 1000) #how often to print updates
    run.mcmc.from.cache(dir=CACHE.DIR,
                        chains=3,
                        update.frequency = 1000) #how often to print updates
    run.mcmc.from.cache(dir=CACHE.DIR,
                        chains=4,
                        update.frequency = 1000) #how often to print updates

    mcmc = assemble.mcmc.from.cache(dir=CACHE.DIR)
}

#If you want the MCMC to cache to the disk as it goes, use arguments
# cache.frequency and cache.dir

##-----------------------------##
##-- DIAGNOSTICS ON THE MCMC --##
##-----------------------------##

print(trace.plot(mcmc)) + ggtitle("Trace Plot")
print("Gelman's Rhat statistic for convergence:")
print(get.rhats(mcmc))

##-----------------------------------------------------------------##
##-- PROFIT! (ie, inference on the posterior set of simulations) --##
##-----------------------------------------------------------------##

simset = extract.simset(mcmc)
#This S4 object now contains the set of all simulations from the MCMC
# We can do fun stuff like this:

#Get the model simulated true notifications
extract.model.notifications.by.year = function(sim)
{
    rv = sim$simulated.notifications
    names(rv) = as.character(YEARS)
    rv
}
notifications.distribution = extract.simset.distribution(simset, extract.model.notifications.by.year)

#This just tells us what our model said every year
print("Simulated Notifications by Year")
print(cbind(mean.notifications=get.means(notifications.distribution),
        t(get.intervals(notifications.distribution))))

#Or we can plot it
df = as.data.frame(cbind(value=get.means(notifications.distribution),
                         t(get.intervals(notifications.distribution))))
df$year = YEARS
df$type = 'Model'
df = rbind(df,
           data.frame(value=NOTIFICATIONS,
                      lower=NA,
                      upper=NA,
                      year=YEARS,
                      type='Truth'))
print(ggplot(df, aes(year, value, ymin=lower, ymax=upper, color=type, fill=type)) +
    geom_ribbon(alpha=0.2) + geom_line(size=1) + geom_point(size=5) +
        ggtitle("Simulated Notifications (mean, 95% CI) vs Truth") + ylab("Notifications") + xlab('Year'))


#We can extract what the notifications reduction is from 2010 to 2020
extract.notifications.reduction = function(sim)
{
    n.2010 = sim$simulated.notifications[1]
    n.2020 = sim$simulated.notifications[11]

    (n.2010-n.2020)/n.2010
}

reduction.distribution = extract.simset.distribution(simset, extract.notifications.reduction)
print("Reduction in Notifications:")
print(c(mean.reduction=get.means(reduction.distribution), get.intervals(reduction.distribution)))
