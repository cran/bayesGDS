rm(list=ls())

library(plyr)
library(Matrix)
library(mvtnorm)
library(trustOptim)
library(bayesGDS)

source("ex_funcs.R")

set.seed(123)


## Set parameters for the sampling algorithm
M <- 20000
ds.scale <- 0.8

n.draws <- 100
max.AR.tries <- 20000


## Simulate data and set priors
N <- 10
k <- 2
T <- 40

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k,0,10)
Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)
X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) 
B <- t(rmvnorm(N, mean=mu, sigma=Omega))
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- laply(log.p, function(q) return(rbinom(1,T,exp(q))))

nvars <- N*k + k
start <- rnorm(nvars) ## random starting values

hess.struct <- get.hess.struct(N, k)  ## sparsity structure of Hessian
 
opt <- trust.optim(start, fn=log.post,
                   gr = get.grad,         
                   hess.struct = hess.struct,
                   method = "SparseFD",
                   control = list(
                     report.freq=1L,         
                     maxit=1000L,
                     function.scale.factor = as.numeric(-1),                           
                     preconditioner=1L
                     ),
                   Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                   )

post.mode <- opt$solution
hess <- opt$hessian


log.c1 <- opt$fval
chol.hess <- t(chol(-hess))

prop.params <- list(mu = post.mode,
                    chol.prec = sqrt(ds.scale)*chol.hess
                    )


log.c2 <- get.log.dens.MVN(post.mode, prop.params)
log.const <- log.c1 - log.c2

cat("Collecting GDS Proposal Draws\n")


draws.m <- as(draw.MVN.proposals(M,prop.params),"matrix")
log.post.m <- aaply(draws.m, 2,log.post,
                    Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma,
                    .parallel=FALSE, .progress="text")
log.prop.m <- get.log.dens.MVN(draws.m,params=prop.params)
log.phi <- log.post.m - log.prop.m +log.c2 - log.c1


cat("Are any log.phi > 0?  ",any(log.phi>0),"\n")

browser()

cat("Generating DS draws - accept-reject phase\n")
draws <- get.GDS.draws(n.draws = n.draws,
                       log.phi=log.phi,
                       log.const = log.const,
                       log.post.func = log.post,
                       draw.prop.func = draw.MVN.proposals,
                       prop.log.dens.func = get.log.dens.MVN,
                       prop.params = prop.params,
                       max.tries=max.AR.tries,
                       max.batch = 20,
                       est.acc.rate=.05,
                       debug=FALSE,
                       report.freq=10,
                       Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                       )

LML <- get.LML(draws$counts, log.phi, log.const)

