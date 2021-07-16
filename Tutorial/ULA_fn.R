library(mvtnorm)

ULA_fn <- function(d,initial,h,cov_ULA,iters,der_loglike,der_logprior,options,varNames){
  
  samples <- d_loglike <- d_logprior <- matrix(NaN,nrow=iters,ncol=d)
  samples[1,] <- initial
  
  temp <- der_loglike(samples[1,],options)
  d_loglike[1,] <- temp$der_loglike
  
  temp <- der_logprior(samples[1,],options)
  d_logprior[1,] <- temp$der_logprior

  for (i in 1:(iters-1)){
    mymean <- samples[i,] + h^2/2*cov_ULA%*%(d_loglike[i,]+d_logprior[i,])
    if (d==1){
      samples[i+1,] <- rnorm(n=1, mymean, h*sqrt(cov_ULA))
    } else {
      samples[i+1,] <-  rmvnorm(n=1, mymean, h^2*cov_ULA) #Multivariate normal random walk for proposal
    }
    
    temp <- tryCatch({
      der_loglike(samples[i+1,],options)
    } , error = function(err){return (NaN)} )
    
    if (any(is.na(temp))){
      d_loglike[i+1,] <- rep(0,d)
      print(sprintf("Error in likelihood evaluation at iteration %d.",i))
    } else{
      d_loglike[i+1,] <- temp$der_loglike
    }
    
    temp <- der_logprior(samples[i+1,],options)
    d_logprior[i+1,] <- temp$der_logprior

  }
  samples <- as.mcmc(samples)
  varnames(samples) <- varNames
  
  return(list(samples=samples,der_loglike=d_loglike,der_logprior=d_logprior))
}
