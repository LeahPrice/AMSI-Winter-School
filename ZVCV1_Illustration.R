samples <- rnorm(100,mean=5,sd=1)
grads <- -(samples-5)
integrand <- samples

X <- cbind(rep(1,100),grads)
betahat <- solve(t(X)%*%X,t(X)%*%integrand)
Vanilla <- mean(integrand)
ZVCV1 <- mean(integrand - X%*%betahat) + betahat[1]
print(Vanilla)
print(ZVCV1)
