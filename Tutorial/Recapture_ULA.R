
args<-commandArgs(TRUE)
print("begin")
print(args)
print("end")
N_index <- eval( parse(text=args[1]) )

print(sprintf("N index is %d.\n",N_index))


mywd <- "/home/users/southl/BSS/final_BetterTuned"
#setwd(mywd)
mysd <- "/storage/users/southl/BSS/ULA_BetterTuned"
mylib.loc <- "/home/users/southl/BSS/pkg/"


source('ULA_fn.R')

load("recapture_ULA_bettertuning.RData")

d <- NROW(cov_ula)

N_set <- c(50,100,250,500,750,1000,2500,10000,1000000)  # Needs to be large enough for k folds
N <- N_set[N_index]

x_all <- list()
u_all <- list()
t_simulate <- rep(NaN,100)
for (k in 1:100){
  set.seed(100*(N_index-1) + k)
  
  #### Getting data
  tic <- proc.time()
  
  myULA <- ULA_fn(d,initial,h,cov_ula,N,der_loglike,der_logprior,options)
  
  t_simulate[k] <- proc.time()[1] - tic[1]
  
  x_all[[k]] <- myULA$samples # samples
  u_all[[k]] <- myULA$der_loglike + myULA$der_logprior # derivatives
  
  print(k)
}

save(x_all, u_all, t_simulate, file = sprintf("%s/Recapture_ULA_%d.RData",mysd,N), envir = environment())
print(sprintf("N=%d is done",N))
