kernal_IMQ <- function(x, beta, band_width){
  s_band_width <- sqrt(band_width)
  x = x / s_band_width
  
  k0 = (1+x)^beta
  k1 <- beta*(1+x)^(beta - 1) / s_band_width
  k2 <- beta*(beta-1)*(1+x)^(beta - 2) / band_width
  
  res <- list(k0 = k0,k1 = k1,k2 = k2)
  return(res)
}

imqKSD <- function(samples, grads, band_width = 1){
  n <- nrow(samples)
  d <- ncol(samples)
  
  norms <- rowSums(samples^2)
  norms <- norms*matrix(1, n, n)
  mu <- tcrossprod(samples, samples)
  norms <- norms - mu
  norms <- norms + t(norms)
  norms <- norms - diag(diag(norms))
  
  sample_grad <- rowSums(samples*grads)
  sample_grad <- sample_grad*matrix(1, n, n)
  mu <- tcrossprod(samples, grads)
  sample_grad <- sample_grad - mu
  sample_grad <- sample_grad + t(sample_grad)
  sample_grad <- sample_grad - diag(diag(sample_grad))
  
  k_res <- kernal_IMQ(norms, -0.5, band_width)
  k_0 <- k_res$k0; k_1 <- k_res$k1; k_2 <- k_res$k2
  k_gram <- -2*d*k_1 - 4*norms*k_2
  k_gram <- k_gram - 2*k_1*sample_grad
  k_gram <- k_gram + tcrossprod(grads, grads)*k_0
  
  res <- sqrt(sum(k_gram))/n
  return(res)
}