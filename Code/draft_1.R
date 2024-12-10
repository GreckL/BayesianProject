### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("bayesmeta")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library(tidyverse)        # Data manipulation 

rm(list=ls())
set.seed(12345)

### 1 - Simulating data ###
#n = 20
#tau_simul = runif(n,min = -1,max = 1)
#sigma_simul = runif(n,min = .5,max = 5)

### Start with frequentist weighted MLE
# Parameters
s=50         # Number of simulations
nstud = 5   # Maximum number of studies

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
                     char = "=")   # Character used to create the bar

for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0,max = 4*theta_true)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul)
      }
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian bagger
    bg[i] <- baggr(df_pooled, "rubin", pooling = "full")
  }
  MLE_mean[n] <- sum(MLE_est)/s
  MLE_MSE[n] <- (MLE_mean[n]-theta_true)^2
  
  # bayesian
  tau_columns <- lapply(bg, function(df) df[["tau"]])
  BAY_MSE[n] <- (sum(unlist(tau_columns))/s-theta_true)^2
  
  setTxtProgressBar(pb, t)
  
}
close(pb)


# Plot the MSEs to compare
df <- data.frame(x=1:nstud, MLE_MSE, BAY_MSE)

df_long <- df %>%
  pivot_longer(cols = c(MLE_MSE, BAY_MSE), names_to = "Series", values_to = "Value")

# Plot with ggplot2
p <- ggplot(df_long, aes(x = x+1, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +  # Add lines
  geom_point(size = 3) +   # Add points
  geom_text(aes(label = round(Value, 1)), vjust = -1, size = 3.5) +  # Add labels
  theme_minimal() +        # Use a minimal theme
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "MSE Comparison",
    x = "Number of Studies",
    y = "MSE"
  )

# Display the plot
x11()
print(p)


