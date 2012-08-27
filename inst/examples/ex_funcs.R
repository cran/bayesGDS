
log.post <- function(pars, Y, X, inv.Omega, inv.Sigma) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=FALSE)
  mu <- pars[(N*k+1):(N*k+k)]

  
  bx <- colSums(X * beta)
  
  log.p <- bx - log1p(exp(bx))
  log.p1 <- -log1p(exp(bx))
  
  data.LL <- sum(lchoose(T,Y) + Y*log.p + (T-Y)*log.p1)
  
  Bmu <- apply(beta, 2, "-", mu)

  normConst <- (N*log(det(inv.Sigma)) + log(det(inv.Omega))- (k*N+1)*log(2*pi))/2
  
  prior <- -0.5 * sum(diag(tcrossprod(Bmu) %*% inv.Sigma))
  hyp <- -0.5 * t(mu) %*% inv.Omega %*% mu
  res <- normConst + data.LL + prior + hyp
  return(as.numeric(res))
  
}

dlog.f.db <- function(pars, Y, X, inv.Omega, inv.Sigma) {

  beta <- matrix(pars[1:(N*k)], k, N, byrow=FALSE)
  mu <- pars[(N*k+1):length(pars)]
  bx <- colSums(X * beta)

  p <- exp(bx)/(1+exp(bx))

  tmp <- Y - T*p

  dLL.db <- apply(X,1,"*",tmp)
    
  Bmu <- apply(beta, 2, "-", mu)
  dprior <- -inv.Sigma %*% Bmu
  
  res <- t(dLL.db) + dprior

  return(as.vector(res))
 
}

dlog.f.dmu <- function(p, Y, X, inv.Omega, inv.Sigma) {

  beta <- matrix(p[1:(N*k)], k, N, byrow=FALSE)
  mu <- p[(N*k+1):length(p)]
  Bmu <- apply(beta, 2, "-", mu)

  res <- inv.Sigma %*% (rowSums(Bmu)) -  inv.Omega %*% mu
  return(res)
}

get.grad <- function(p, Y, X, inv.Omega, inv.Sigma) {

  q1 <- dlog.f.db(p, Y, X, inv.Omega, inv.Sigma)
  q2 <- dlog.f.dmu(p, Y, X, inv.Omega, inv.Sigma)
  res <- c(q1, q2)
  return(res)
}


get.hess.struct <- function(N, k) {

  B1 <- kronecker(Diagonal(N),Matrix(TRUE,k,k))
  B2 <- Matrix(TRUE,k,N*k)
  B3 <- Matrix(TRUE,k,k)
  H <- cBind(rBind(B1,B2),rBind(t(B2),B3))
  res <- Matrix.to.Coord(H)
  return(res)
  
}
