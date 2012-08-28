\name{get.GDS.draws}
\alias{get.GDS.draws}
\title{Collect draws from the target posterior distribution.}
\description{
Runs the accept-reject phase of the Braun and Damien (2012) Generalized Direct Sampling algorithm.
	}
\usage{
get.GDS.draws(n.draws, log.phi, log.const, log.post.func, draw.prop.func,	
 prop.log.dens.func, prop.params, max.tries=50000, max.batch=1000, est.acc.rate=.1, ...,
 debug=FALSE, report.freq=1)	
}

\arguments{
\item{n.draws}{number of draws to take from the target posterior density.}
\item{log.phi}{Vector of log.phi, as computed from the proposal draws.} 
\item{log.const}{Difference between the log posterior density and the log proposal density, evaluated at the posterior mode.  Using the terminology in Braun and Damien (2012), this is log(c_1) - log(c_2).}
\item{log.post.func}{Function that returns the log posterior density.  Function should take the parameter vector as the first argument.}
\item{draw.prop.func}{Function to draw from the proposal distribution.}
\item{prop.log.dens.func}{Function to return log density of the proposal distribution.  Function should take the parameter vector as the first argument.}
\item{prop.params}{Object (list or vector) to be passed to both draw.prop.func and prop.log.dens.func.  Contains parameters for the proposal distribution.  See details.}
\item{max.tries}{Maximum number of proposal draws to try, without a success.  This prevents the routine from being stuck in an endless loop.}
\item{max.batch}{The algorithm draws proposals in batches, instead of one-by-one.  This argument controls the maximum batch size.  It is used to manage memory consumption.}
\item{est.acc.rate}{Estimated acceptance rate.  This is used to determine the size of the initial batch of proposal draws.  The reason we do this, instead of drawing proposals one by one, is to take advantage of vectorization in the BLAS.}
\item{...}{Additional parameters to be passed to log.post.func.}
\item{debug}{If TRUE, will print various tracing messages to monitor progress of the algorithm.}
\item{report.freq}{Frequency at which to report the sample number that is being currently drawn.  Will be between 1 and n.draws.  Used to monitor progress of the algorithm.}
}
\value{
a list with the following elements:
\item{draws}{A matrix with each draw in a row, and each parameter in a column}
\item{counts}{The number of attempts that it took to get an accepted draw.  The accepted draw counts, so the count will always be at least 1.}
\item{gt.1}{A vector that indicates if the phi for that draw was greater than 1.  Available as a diagnostic.  Normally, these should all be FALSE.  Any values of TRUE suggest that a change in proposal density might be warranted.}
\item{log.post.dens}{A numeric vector.  Log posterior density for each draw.}
\item{log.prop.dens}{A numeric vector. Log of the proposal density for each draw.}
\item{log.thresholds}{Vector of threshold draws (log u) from the accept-reject algorithm.  Sorted in ascending order.}
\item{log.phi}{A numeric vector.  Value of log.phi for the accepted draws.}
}
\details{
For use in conjunction with the Braun and Damien (2012) Generalized Direct Sampling algorithm. 
The prop.params argument contains parameters that are passed to the functions that sample from, and compute the log density of, the proposal distribution. If you are using the draw.MVN.proposals and get.log.dens.MVN functions that are included in this package, there are two elements in this list.  The element mu is the mean of the MVN, and chol.prec is the Cholesky decomposition of the precision matrix (inverse of the covariance matrix).  
}






 


