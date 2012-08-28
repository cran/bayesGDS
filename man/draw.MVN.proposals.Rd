\name{draw.MVN.proposals}
\alias{draw.MVN.proposals}
\title{Draw from MVN, given mean and the Cholesky decomposition of the precision matrix}
\description{
Returns a single multivariate normal draw.  Parameterization is a mean vector, and the Cholesky decomposition of the precision matrix.
	}
\usage{
draw.MVN.proposals(n, params)
}
\arguments{
	\item{n}{number of draws}
	\item{params}{a list with the following elements:
\describe{
\item{mu}{Mean}
\item{chol.prec}{A lower triangular matrix that contains the Cholesky decomposition of the precision matrix.}
}
}
}
\details{
The precision matrix is simply the inverse of the covariance matrix.  Sometimes it is more efficient to work with precision matrices.  This function allows you to compute the Cholesky decomposition only once, and then call this function multiple times to get multiple MVN draws.  For stability and efficiency, we recommend passing chol.prec as an object of class triangularMatrix from the Matrix package.  
}
\value{
a matrix with n columns, with each draw in a column. The number of rows is equal to the length of params$mu.
}







 


