### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("bayesmeta")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library(tidyverse)        # Data manipulation 
library(rstan)
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
BAY_MSE <- c()
bg <- c()
BAY_tau <- c()
bg_mse <- c()
MSE_hat <- 0
MSE_BAY <- c()
theta_true <- 1

# keep track
pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = s, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s)
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0,max = 4*theta_true)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul[t])
      }
   
    # Frequentist
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    ## Using partial
    # Bayesian bagger
    bg <- baggr(df_pooled, model="rubin", pooling = "partial")
    # Extract the Stan fit object
    fit_object <- bg$fit
    # Get the summary table from the Stan fit object
    fit_summary <- summary(fit_object)$summary
    # Convert it to a data frame for easier manipulation
    fit_summary_df <- as.data.frame(fit_summary)
    # Get row names containing "theta_k"
    theta_k_rows <- fit_summary_df[grep("theta_k", rownames(fit_summary)), ]
    theta_SE <- (theta_k_rows$mean-1)^2
    BAY_MSE <- sum(unlist(theta_SE))/n
    #if (s == 1) {
    #  MSE_hat = BAY_MSE
    #}
    #else{
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s)
    #}
    #*BAY_Nstud_MSE[i] <- sum(BAY_MSE)/n
    # Print the theta_k rows
    #print(theta_k_rows)
    
  }
  MLE_mean[n] <- sum(MLE_est)/s
  MLE_MSE[n] <- (MLE_mean[n]-theta_true)^2
  
  # bayesian
  MSE_BAY[n] = MSE_hat
  #tau_columns <- lapply(bg, function(df) df[["tau"]])
  #BAY <- ba
  #BAY_MSE[n] <- ((sum(unlist(tau_columns))/n)/s-theta_true)^2
  #BAY_tau[i] <- hypermean(bg[i], transform = NULL,
  #                        interval = 0.95,
  #                        message = FALSE,
  #                        summary = TRUE)
  
  #setTxtProgressBar(pb, t)
  #bg_mse$n = (tau_columns[[n]]-1)^2
  #BAY_MSE[n] = sum(unlist(bg_mse))/n
  
}
close(pb)


# Plot the MSEs to compare
df <- data.frame(x=1:nstud, MLE_MSE, MSE_BAY)

df_long <- df %>%
  pivot_longer(cols = c(MLE_MSE, MSE_BAY), names_to = "Series", values_to = "Value")

# Plot with ggplot2
p <- ggplot(df_long, aes(x = x, y = Value, color = Series, group = Series)) +
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


