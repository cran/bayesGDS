\name{get.log.dens.MVN}
\alias{get.log.dens.MVN}
\title{Evalaute a multivariate normal log density,  given mean and the Cholesky decomposition of the precision matrix}
\description{
Returns the log density of a multivariate normal density, evaluated at x.  Parameterization is a mean vector, and the Cholesky decomposition of the precision matrix.
	}
\usage{
get.log.dens.MVN(x, params)
}

\arguments{
	\item{x}{vector at which to evaluate the density}
	\item{params}{a list with the following elements:
\describe{
	\item{mu}{Mean}
	\item{chol.prec}{A lower triangular matrix that contains the Cholesky decomposition of the precision matrix.}
}
}
}
\value{
a numeric.  the log density of an MVN, evaluated at x.
}
\details{
The precision matrix is simply the inverse of the covariance matrix.  Sometimes it is more efficient to work with precision matrices.  This function allows you to compute the Cholesky decomposition only once, and then call this function multiple times to compute densities for multiple x's.
}






 


