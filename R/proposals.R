
draw.MVN.proposals <- function(n, params) {

  ## n = number of draws to collect.  Each column is a draw
  ## params is a list:
  ##   mu - mean of the MVN
  ##   chol.prec - lower-triangular Cholesky of the precision matrix
  
  k <- length(params$mu)
  x <- matrix(rnorm(n*k),k,n)  
  y <- as(solve(Matrix:::t(params$chol.prec),x),"matrix")

  res <- y+matrix(params$mu,k,n,byrow=FALSE)

  return(res)

}

get.log.dens.MVN <- function(x, params) {

 
  ## INPUT:
  ##    x : point at which density is evaluated
  ##    params:  a list with the following elements:
  ##        mu : the mean of the MVN
  ##        chol.prec : a lower triangular Matrix object, containing
  ##                      the Cholesky factor of the PRECISION matrix

  ## OUTPUT:
  ##    log density of MVN evaluated at x
  
  k <- NROW(x)
  M <- NCOL(x)
  
  xmu <- x-matrix(params$mu,nrow=k, ncol=M,byrow=FALSE)

  ## something weird is going on with namespaces
  y <- t(params$chol.prec) %*% xmu

  log.dens <- -k*log(2*pi)/2 + sum(log(diag(params$chol.prec))) - colSums(y*y)/2
 
  return(as.numeric(log.dens))
}

