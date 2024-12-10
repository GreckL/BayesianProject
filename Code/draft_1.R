### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("bayesmeta")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 


rm(list=ls())
set.seed(12345)

### 1 - Simulating data ###
#n = 20
#tau_simul = runif(n,min = -1,max = 1)
#sigma_simul = runif(n,min = .5,max = 5)

### Start with frequentist weighted MLE
s=100
#MLE_est <-rep(NA,s)

theta_simul <- c()
sigma_simul <- c()
MLE_est <- c()
MLE_mean <- c()
MLE_MSE <- c()
bg <- c()
theta_true <- 1

# keep track
pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = s, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for (n in 1:100) {      # Number of studies
  for (i in 1:s) {      # Number of simulations
    for (t in 2:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0,max = 4*theta_true)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul)
      }
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian bagger
    bg[i] <- baggr(df_pooled, pooling = "none")
  }
  MLE_mean[n] <- sum(MLE_est)/s
  MLE_MSE[n] <- (MLE_mean[n]-theta_true)^2
  
  # bayesian
  
  setTxtProgressBar(pb, t)
  
}
close(pb)



ggplot(MLE_MSE)
plot( MLE_MSE)

df_pooled <- data.frame("tau" = tau_simul,
                        "se"  = sigma_simul)
bg <- baggr(df_pooled, pooling = "full")




