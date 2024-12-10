### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("bayesmeta")
library("baggr")
library("bayesplot")

rm(list=ls())
set.seed(12345)

### 1 - Simulating data ###
n = 20
tau_simul = runif(n,min = -1,max = 1)
sigma_simul = runif(n,min = .5,max = 5)


### Start with frequentist weighted MLE
s=10000
MLE_est <-rep(NA,s)

for i in 1:s {
  tau_simul = runif(n,min = -1,max = 1)
  sigma_simul = runif(n,min = .5,max = 5)
  for t in 1:n {
    MLE_est[i]<-sum(tau_simul)
  }
}


df_pooled <- data.frame("tau" = tau_simul,
                        "se"  = sigma_simul)
bg <- baggr(df_pooled, pooling = "full")




