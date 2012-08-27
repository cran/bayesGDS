
get.LML <- function(counts, log.phi, log.const, bits=128) {

  if ( any(!is.finite(counts)) | any(counts<=0) ) {
    stop ("Error in get.LML:  all counts must be finite and positive")
  }

  if ( any(!is.finite(log.phi)) | any(log.phi>0) ) {
    stop ("Error in get.LML:  all values of log.phi must be non-positive")
  }

  
  ar <- 1/mean(counts)
  M <- length(log.phi)
  ord <- order(log.phi,decreasing=TRUE)
  v <- mpfr(log.phi[ord], bits)
  max.v <- max(v)
  z <- sum((2*seq(1,M)-1) * exp(v-max.v))

  LL <- max.v + log.const - log(ar) - 2*log(M) + log(z)
  return(as.numeric(LL))
 

}

