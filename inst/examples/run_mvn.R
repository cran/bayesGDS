
## Section A: Preamble

rm(list=ls())
gc()
library(bayesGDS)
detach(package:bayesGDS,unload=TRUE)
library(Rmpfr)
library(Matrix)
library(foreach)
library(doParallel)
library(plyr)
library(trustOptim)
library(bayesGDS)

run.par <- TRUE
if(run.par) registerDoParallel(cores=12) else registerDoParallel(cores=1)

set.seed(12345)

## Section B:  Functions for LL, gradient and Hessian structure


get.log.post.HMVN <- function(x, data) {
  f <- .Call("HierMVN_get_f",as.vector(x),data) ## C function called from trustOptim package
  return(f)

}

get.log.post.HMVN.grad <- function(x, data) {
  df <- .Call("HierMVN_get_fdf", as.vector(x), data)$grad
  return(df)
}



get.hess.struct.HMVN <- function(nvars, n.obs, k){

  ## construct hessian structure -- block-arrow structure - lower.tri only

  b.range <- 1:(k*n.obs)
  
  HS <- Matrix(0,nrow=nvars, ncol=nvars)
  
  HS[1:(k*n.obs),1:(k*n.obs)] <- kronecker(Diagonal(n.obs),Matrix(1,k,k)) ## check conditional dependence
  HS[(k*n.obs+1):nvars,] <- 1

  HS <- as(tril(HS),"TsparseMatrix")
  iRow <- as.integer(HS@i+1)
  jCol <- as.integer(HS@j+1)

  return(list(iRow=iRow, jCol=jCol))
}

## Section C:  Set parameters for the sampling algorithm
n.draws <- 20
M <- 5000
max.AR.tries <- 2500
ds.scale <- 0.71
get.log.prop <- get.log.dens.MVN
draw.prop <- draw.MVN.proposals

## Section D:  Simulate data
n.obs <- 100
T <- 6
k <- 3


X <- rbind(rep(1,n.obs),matrix(rnorm((k-1)*n.obs),nrow=k-1,ncol=n.obs))
beta <- seq(-5,5,length=k)
sig <- 1
mu <- beta %*% X
Y <- t(laply(mu,function(x) return(rnorm(T,mean=x,sd=sig))))
data.list <- list(X=X, Y=Y)

## Section E:  Set priors

nvars <- as.integer(n.obs*k + k + (k*(k+1)/2))
mu.prior.mean <- rep(0,k)
mu.prior.chol.prec <- t(chol(5*diag(k)))
nu <- as.numeric(k+6)
G.mean <- 2*diag(k)
chol.inv.mean.G <- t(chol(solve(G.mean)))

prior.list <- list(mu.prior.mean=mu.prior.mean,
                   mu.prior.chol.prec=mu.prior.chol.prec,
                   nu=nu,
                   chol.inv.mean.G=chol.inv.mean.G
                   )

## Section F:  Find posterior mode

params.list <- list(data=data.list, priors=prior.list)
hs.list <- get.hess.struct.HMVN(nvars, n.obs, k)  ## indexing starts at 1


B.start <- rnorm(n.obs*k)
mu.start <- mu.prior.mean
G.start <- vech(chol.inv.mean.G)

startX <- c(B.start, mu.start, G.start)
names(startX) <- paste("x_",1:nvars,sep="")

control.list <- list(step.size=as.numeric(20),
                     tol=sqrt(.Machine$double.eps),
                     prec=sqrt(.Machine$double.eps),
                     report.freq=as.integer(10),
                     report.level=as.integer(5),
                     report.precision=as.integer(5),
                     maxit=as.integer(4000),
                     contract.factor = 0.5,
                     expand.factor = 3,
                     contract.threshold = 0.15,
                     expand.threshold.ap = 0.85,
                     expand.threshold.radius = 0.85,
                     function.scale.factor = as.numeric(-1),
                     hessian.refresh.freq=as.integer(1L),
                     precond.refresh.freq=as.integer(100L),
                     quasi.newton.method=as.integer(0),  ## 1 = SR1, 2 = BFGS
                     preconditioner=as.integer(1),  ## 0=none, 1=diagonal 2=ssor 3=cholesky
                     fd.method=as.integer(1)
                     )


opt <- trust.optim(startX, fn=get.log.post.HMVN, gr=get.log.post.HMVN.grad,
                   method="SparseFD",
                   hess.struct=hs.list,
                   control=control.list,
                   data=params.list)

post.mode <- opt$solution
hess <- opt$hessian
chol.hess <- t(chol(-hess))


## Section G:  The Generalized Direct Sampler

get.f <- function(x,...) return(get.log.post.HMVN(x,params.list))

prop.params <- list(mu = post.mode,
                    chol.prec = sqrt(ds.scale)*chol.hess
                    )

log.c1 <- opt$fval
log.c2 <- get.log.prop(post.mode, prop.params)
log.const <- log.c1 - log.c2

cat("Collecting GDS Proposal Draws\n")


draws.m <- as(draw.MVN.proposals(M,prop.params),"matrix")
log.post.m <- aaply(draws.m, 2,get.f, .parallel=TRUE, .progress="none")
log.prop.m <- get.log.dens.MVN(draws.m,params=prop.params)
log.phi <- log.post.m - log.prop.m +log.c2 - log.c1


cat("Are any log.phi > 0?  ",any(log.phi>0),"\n")

cat("Generating DS draws - accept-reject phase\n")
draws <- get.GDS.draws(
                       n.draws = n.draws,
                       log.phi=log.phi,
                       log.const = log.const,
                       log.post.func = get.f,
                       draw.prop.func = draw.MVN.proposals,
                       prop.log.dens.func = get.log.dens.MVN,
                       prop.params = prop.params,
                       max.tries=max.AR.tries,
                       max.batch = 20,
                       est.acc.rate=.05,
                       debug=FALSE,
                       report.freq=1)

LML <- get.LML(draws$counts, log.phi, log.const)



