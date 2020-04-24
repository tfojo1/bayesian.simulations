
library(mcmc.sim)

if (!grepl('jheem', getwd()))
    setwd("../../../Ending HIV/jheem/")
source('../code/calibration/run_for_parameters.R')

if (1==2)
{
    components = setup.components.for.msa(BALTIMORE.MSA, population.year=2007)

    starting.values = c(hiv.mortality=0.0869*.8,
                        peak.trate.mult=1.8,
                        idu.trate=0.0025,
                        sexual.trate=0.01*.95)

    starting.values = c(hiv.mortality=0.0849, peak.trate.mult=2.4942, idu.trate=0.0018, sexual.trate=0.008)

    sim.function = function(pp)
    {
        parameters=c(hiv.mortality=0.0869*.8,
                     black.black.oe=3.75,
                     hispanic.hispanic.oe=2.2,
                     other.other.oe=1.55,
                     black.susceptibility.rr=4.5,
                     hispanic.susceptibility.rr=2.5,
                     msm.transmission=15,
                     peak.trate.mult=1.8,
                     idu.trate=0.0025,
                     sexual.trate=0.01*.95)

        parameters[names(pp)] = pp

        sim = run.jheem.for.parameters(components=components, parameters=parameters)

        print(plot.single.calibration.total(sim))

        sim
    }

    log.prior = function(params){0}
    likelihood.total = create.joint.likelihood.function(create.likelihood.function(data.type='new',by.total=T,
                                                                                   year.to.year.correlation=0.5),
                                                        create.likelihood.function(data.type='prevalence',by.total=T,
                                                                                   year.to.year.correlation=0.5))

    ctrl = create.adaptive.blockwise.metropolis.control(var.names=names(starting.values),
                                                       simulation.function = sim.function,
                                                       log.prior.distribution = log.prior,
                                                       log.likelihood = likelihood.total,
                                                       initial.covariance.mat = diag(starting.values/40)^2)

    mcmc = run.initial.chain(ctrl, start.parameters = starting.values,
                                         n.iter=200, update.frequency = 1)



    ctrl.plain = create.metropolis.control(var.names=names(starting.values),
                                           simulation.function = sim.function,
                                           log.prior.distribution = log.prior,
                                           log.likelihood = likelihood.total,
                                           proposal.sds = starting.values/40,
                                           use.adaptive.scaling = T)


    mcmc = run.initial.chain(ctrl.plain, start.parameters = starting.values,
                             n.iter=200, update.frequency = 1)

    ctrl = create.adaptive.blockwise.metropolis.control(var.names=names(starting.values),
                                                        simulation.function = sim.function,
                                                        log.prior.distribution = log.prior,
                                                        log.likelihood = likelihood.total,
                                                        initial.covariance.mat = diag(starting.values/40)^2,
                                                        use.adaptive.covariance = F)

    mcmc = run.initial.chain(ctrl, start.parameters = starting.values,
                             n.iter=200, update.frequency = 1)
}

#really good params:
#hiv.mortality=0.0849, peak.trate.mult=2.4942, idu.trate=0.0018, sexual.trate=0.008
