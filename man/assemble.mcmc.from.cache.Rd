% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_interface.R
\name{assemble.mcmc.from.cache}
\alias{assemble.mcmc.from.cache}
\title{Assemble an MCMC object from iterations stored in a cache}
\usage{
assemble.mcmc.from.cache(
  dir,
  allow.incomplete = F,
  chains = NULL,
  remove.cache.when.done = F
)
}
\arguments{
\item{dir}{The directory where the cache is stored}

\item{allow.incomplete}{If FALSE, throws an error if the cache has not finished running. If TRUE, assembles an MCMC object out of all chains where any iterations have been run past the burn period, with a number of iterations equal to the mininum number of iterations run in any chain past the burn period}

\item{chains}{The numbers of the chains to run. If passed NULL (the default), runs all chains}
}
\value{
An \link{mcmc.sim} object
}
\description{
Assemble an MCMC object from iterations stored in a cache
}
