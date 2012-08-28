


draw.thresholds <- function(log.phi, n.draws, bits=128) {

  ## INPUTS:
  ##   log.phi - log of \Phi (should all be negative)
  ##   n.draws - number of draws required from the posterior
 

  ## OUTPUTS:
  ##   log.u - threshold draws, using the log transformation in Braun and Damien, for use in GDS A-R step

  if (any(log.phi>0) ) stop("all values of log phi must be non-positive")

  M <- length(log.phi)
  v <- mpfr(c(-log.phi[order(log.phi,decreasing=TRUE)],Inf), precBits=bits)
  v1 <- min(v)

  p.int <- as.numeric(exp(log(seq(1,M)) -v1 + log(exp(v1-v[1:M]) - exp(v1-v[2:(M+1)]))))

  j <- sample(1:M,n.draws,replace=TRUE,prob=p.int)

  u <- runif(n.draws)
  d <- as.numeric(v[j] - log(1-u + u*exp(-v[j+1]))) ## draw from right-censored exponential

  return(d[order(d,decreasing=TRUE)])
}




get.GDS.draws<- function(n.draws, log.phi, log.const, log.post.func, draw.prop.func,
                         prop.log.dens.func, prop.params, max.tries=50000, max.batch=1000,
                         est.acc.rate=.1,..., debug=FALSE, report.freq=1) {
  

## serial, batched version.  Idea is to run this whole function in parallel, and collect draws at the end
  
  v <- draw.thresholds(log.phi, n.draws)
  nvars <- length(prop.params$mu)


  draw.mat <- matrix(nrow=nvars, ncol=n.draws)
  counts <- vector("integer",length=n.draws)
  gt.1 <- vector("integer",length=n.draws)
  log.post.dens <- vector("numeric",length=n.draws)
  log.prop.dens <- vector("numeric",length=n.draws)

  if (debug) cat("-Drawing batch of proposals\n")
  batch.size <- min(max(ceiling(n.draws/est.acc.rate),10),max.batch)
  props <- as(draw.prop.func(batch.size, prop.params),"matrix")
  if (debug) cat("Computing proposal densities\n")
  log.prop.x <- prop.log.dens.func(props, prop.params)

  if (debug) cat("Computing log posterior at all proposals\n")
  log.post.x <- rep(0,batch.size)

  for (j in 1:batch.size) {
    log.post.x[j] <- log.post.func(props[,j],...)
  }

  crit <- log.post.x - log.prop.x - log.const

  idx <- 1
 
  for (i in 1:n.draws) {
    if ( (i %% report.freq) == 0) cat(" -- draw ",i,"\n")
    stop.loop <- FALSE
    count <- 1
    while (!stop.loop) {
     if (debug) {
       if ((idx %% 500) == 0) cat("\t\tidx = ",idx,"\n")
     }
     
     if (is.finite(crit[idx])) {
       if (-v[i] < crit[idx]) {
         
         draw.mat[,i] <- props[,idx]
         counts[i] <- count
         log.post.dens[i] <- log.post.x[idx]
         log.prop.dens[i] <- log.prop.x[idx]
         stop.loop <- TRUE
         if (crit[idx]>0) gt.1[i] <- TRUE else gt.1[i] <- FALSE
       } else {
         count <- count + 1
         if (count >= max.tries) {
           draw.mat[,i] <- rep(NA,nvars)
           counts[i] <- NA
           log.post.dens[i] <- NA
           log.prop.dens[i] <- NA
           stop.loop <- TRUE
           gt.1[i] <- NA
         }
       }
      } else {
        count <- count + 1
      }
      idx <- idx+1
      if (idx > ncol(props)) {
 
     if (debug) cat("--Drawing a new batch of proposals after draw ",i," and proposal ",idx,"\n")
  
        remaining.draws <- n.draws-i+1
        act.acc.rate <- 1/mean(counts[1:i-1],na.rm=TRUE)
        if (!is.finite(act.acc.rate)) act.acc.rate <- .001
        batch.size <- min(max(ceiling(remaining.draws/act.acc.rate),10),max.batch)
    
     if (debug) cat("---- drawing proposals\n")
        rm(props)
        gc()
        props <- as(draw.prop.func(batch.size, prop.params),"matrix")
     if (debug) cat("Computing proposal densities\n")
        log.prop.x <- prop.log.dens.func(props, prop.params)        


     if (debug) cat("Computing log posterior at all proposals\n")
        log.post.x <- rep(0,batch.size)
        for (j in 1:batch.size) {
          log.post.x[j] <- log.post.func(props[,j],...)
        }

        crit <- log.post.x - log.prop.x - log.const
        idx <- 1
      }      
    } ## end while loop
  }  ## end loop on draws

  return(list(draws=draw.mat,
              counts=counts,
              gt.1=gt.1,
              log.post.dens=log.post.dens,
              log.prop.dens=log.prop.dens,
              log.thresholds=v,
              log.phi=log.post.dens-log.prop.dens-log.const)
         )

}

