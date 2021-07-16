

MALA_fn <- function(d,initial,h,covmala,iters,der_loglike,der_logprior,options,varNames){
  
  samples <- d_loglike <- d_logprior <- matrix(NaN,nrow=iters,ncol=d)
  loglike <- logprior <- rep(NaN,iters)
  samples[1,] <- initial
  
  temp <- der_loglike(samples[1,],options)
  loglike[1] <- temp$loglike
  d_loglike[1,] <- temp$der_loglike
  
  temp <- der_logprior(samples[1,],options)
  logprior[1] <- temp$logprior
  d_logprior[1,] <- temp$der_logprior
  
  for (i in 1:(iters-1)){
    mymean <- samples[i,] + h^2/2*covmala%*%(d_loglike[i,]+d_logprior[i,])
    if (d==1){
      samples_prop <- rnorm(n=1, mymean, h*sqrt(covmala))
    } else {
      samples_prop <-  rmvnorm(n=1, mymean, h^2*covmala) #Multivariate normal random walk for proposal
    }
    
    temp <- tryCatch({
      der_loglike(samples_prop,options)
    } , error = function(err){return (NaN)} )
    
    if (any(is.na(temp))){
      loglike_prop <- -Inf
      d_loglike_prop <- 0
      print("Error in likelihood evaluation")
    } else{
      loglike_prop <- temp$loglike
      d_loglike_prop <- temp$der_loglike
    }
    
    
    temp <- der_logprior(samples_prop,options)
    logprior_prop <- temp$logprior
    d_logprior_prop <- temp$der_logprior
    
    mymean_prop <- t(samples_prop) + h^2/2*covmala%*%(d_loglike_prop+d_logprior_prop)
    
    if (d==1){
      transition_totheta <- dnorm(t(samples[i,]),mymean_prop,h*sqrt(covmala),log=TRUE);
      transition_toprop <- dnorm(samples_prop,mymean,h*sqrt(covmala),log=TRUE);
    }	else{
      transition_totheta <- dmvnorm(t(samples[i,]),mymean_prop,h^2*covmala,log=TRUE);
      transition_toprop <- dmvnorm(samples_prop,mymean,h^2*covmala,log=TRUE);
    }
    
    log_mh <- (loglike_prop - loglike[i]) + logprior_prop - logprior[i] + transition_totheta - transition_toprop
    
    # determine whether to accept or reject
    if (exp(log_mh) > runif(1,0,1)){
      # then accept the proposal
      samples[i+1,] <- samples_prop
      loglike[i+1] <- loglike_prop
      logprior[i+1] <- logprior_prop
      d_loglike[i+1,] <- d_loglike_prop
      d_logprior[i+1,] <- d_logprior_prop
    } else{
      samples[i+1,] <- samples[i,]
      loglike[i+1] <- loglike[i]
      logprior[i+1] <- logprior[i]
      d_loglike[i+1,] <- d_loglike[i,]
      d_logprior[i+1,] <- d_logprior[i,]
    }
    # if (i%%100){
    #   print(i)
    # }
  }  
  samples <- as.mcmc(samples)
  varnames(samples) <- varNames
  
  return(list(samples=samples,loglike=loglike,logprior=logprior,der_loglike=d_loglike,der_logprior=d_logprior))
  
}
recapture_simprior <- function(N,options){
  part_vals <- matrix(runif(N*11),nrow=N,ncol=11) #draw from prior
  part_vals <- log(part_vals/(1-part_vals)) #transform to real line
  return(part_vals)
}

RW_fn <- function(d,initial,h,covmala,iters,der_loglike,der_logprior,options,varNames){
  
  samples <- d_loglike <- d_logprior <- matrix(NaN,nrow=iters,ncol=d)
  loglike <- logprior <- rep(NaN,iters)
  samples[1,] <- initial
  
  temp <- der_loglike(samples[1,],options)
  loglike[1] <- temp$loglike
  d_loglike[1,] <- temp$der_loglike
  
  temp <- der_logprior(samples[1,],options)
  logprior[1] <- temp$logprior
  d_logprior[1,] <- temp$der_logprior
  
  for (i in 1:(iters-1)){
    if (d==1){
      samples_prop <- rnorm(n=1, samples[i,], h*sqrt(covmala))
    } else {
      samples_prop <-  rmvnorm(n=1, samples[i,], h^2*covmala) #Multivariate normal random walk for proposal
    }
    
    temp <- tryCatch({
      der_loglike(samples_prop,options)
    } , error = function(err){return (NaN)} )
    
    if (any(is.na(temp))){
      loglike_prop <- -Inf
      d_loglike_prop <- 0
      print("Error in likelihood evaluation")
    } else{
      loglike_prop <- temp$loglike
      d_loglike_prop <- temp$der_loglike
    }
    
    
    temp <- der_logprior(samples_prop,options)
    logprior_prop <- temp$logprior
    d_logprior_prop <- temp$der_logprior
        
    log_mh <- (loglike_prop - loglike[i]) + logprior_prop - logprior[i]
    
    # determine whether to accept or reject
    if (exp(log_mh) > runif(1,0,1)){
      # then accept the proposal
      samples[i+1,] <- samples_prop
      loglike[i+1] <- loglike_prop
      logprior[i+1] <- logprior_prop
      d_loglike[i+1,] <- d_loglike_prop
      d_logprior[i+1,] <- d_logprior_prop
    } else{
      samples[i+1,] <- samples[i,]
      loglike[i+1] <- loglike[i]
      logprior[i+1] <- logprior[i]
      d_loglike[i+1,] <- d_loglike[i,]
      d_logprior[i+1,] <- d_logprior[i,]
    }
    if (i%%10000==0){
      print(i)
    }
  }  
  samples <- as.mcmc(samples)
  varnames(samples) <- varNames
  
  return(list(samples=samples,loglike=loglike,logprior=logprior,der_loglike=d_loglike,der_logprior=d_logprior))
  
}
