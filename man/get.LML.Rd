\name{get.LML}
\alias{get.LML}
\title{Estimate log marginal likelihood of model}
\description{
Estimate log marginal likelihood of model from GDS output
	}
\usage{
get.LML(counts, log.phi, log.const, bits=128)
}

\arguments{
	\item{counts}{vector of counts of the number of proposals that were generated before accepting a draw.  Length of vector is equal to the number of draws from the posterior.  If the first proposal for a particular posterior draw is accepted, that count is a 1.}
	\item{log.phi}{Numeric vector of draws of log.phi from the proposal draws.}
	\item{log.const}{log_c1 - log_c2}
	\item{bits}{Number of precision bits, as used in the Rmpfr package.}
}
\value{
the estimate log marginal likelihood of the model.
}






 


