\name{draw.thresholds}
\alias{draw.thresholds}
\title{Draw thresholds for the accept-reject stage of the GDS sampling algorithm.}
\description{
Returns a vector of log(u), where u is the threshold to determine if a proposal draw should be accepted as a draw from the target posterior distribution.
	}
\usage{
draw.thresholds(log.phi, n.draws, bits=128)
}

\arguments{
	\item{log.phi}{Vector of log.phi from the proposal draws.  All must be non-positive.}
	\item{n.draws}{an integer.  number of draws to be taken from the target posterior.}
	\item{bits}{Number of precision bits, as used in the Rmpfr package.}
}
\value{
a numeric vector for v = -log.u (the thresholds for the accept-reject stage).
}
\details{
For use in conjunction with the Braun and Damien (2011) Generalized Direct Sampling algorithm.  This is usually not called directly, since it is called from the get.GDS.draws function.
}






 


