### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library("tidyverse")        # Data manipulation 
library("rstan")
library("compiler")

rm(list=ls())
set.seed(12345)

########## Situation 1: DGM s.t. Fixed Effects and no Covariates ###############

# Parameters
s=500         # Number of simulations
nstud = 50   # Maximum number of studies
theta_simul <- c()  #Simulate data set
sigma_simul <- c()  #Simulate data set
MLE_est <- c()    #Vector of simulated Frequentist weighted-pool Estimators
MLE_mean <- c()   #Monte Carlo average 
MLE_MSE <- c()    #Mean Squared Error of Frequentist approach
BAY_MSE <- c()    #Mean Squared Error of Bayesian approach
bg <- c()     #Bayesian Model
BAY_tau <- c()  
bg_mse <- c()   #Bayesian Mean Squared Error
MSE_hat <- 0
MSE_BAY <- c()
theta_true <- 1   #True Fixed Treatment Effect

custom_prior  <- list(hypermean = normal(0,8^2), hypersd = uniform(0,10*0.25))

### Start with partial pooling
for (n in seq(2, nstud, by = 1))  {   
  for (i in 1:s) {      
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = rexp(1, rate = theta_true*4)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul[t])
      }
    # Frequentist
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    ## Using partial pooling from Bayesian bagger
    bg <- baggr(df_pooled, model="rubin", pooling = "partial", prior = custom_prior)
    # Extract the Stan fit object
    fit_object <- bg$fit
    # Get the summary table from the Stan fit object
    fit_summary <- summary(fit_object)$summary
    # Convert it to a data frame for easier manipulation
    fit_summary_df <- as.data.frame(fit_summary)
    # Get row names containing "theta_k"
    theta_k_rows <- fit_summary_df[grep("theta_k", rownames(fit_summary)), ]
    theta_MSE <- (theta_k_rows$mean-1)^2
    BAY_MSE <- sum(unlist(theta_MSE))/n
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s)
    #*BAY_Nstud_MSE[i] <- sum(BAY_MSE)/n
    print(i)
  }
  MLE_mean[n] <- sum(MLE_est)/s
  MLE_MSE[n] <- (MLE_mean[n]-theta_true)^2
  
  # Bayesian
  MSE_BAY[n] = MSE_hat
}

# Plot the MSEs to compare
df <- data.frame(x=1:nstud, MLE_MSE, MSE_BAY)
dflong <- df %>%
  pivot_longer(cols = c(MLE_MSE, MSE_BAY), names_to = "Series", values_to = "Value")

# Plot with ggplot2
p <- ggplot(dflong, aes(x = x, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +  # Thicker lines for clarity
  geom_point(size = 2) +   # Distinct points for readability
  geom_text(aes(label = round(Value, 1)), 
            vjust = -0.7, hjust = 0.7, 
            size = 3.5, family = "sans", color = "black",
            check_overlap = TRUE) +  # Prevent text overlap
  theme_minimal(base_family = "sans") +  # Professional theme
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),  # Lighter gridlines
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("MLE_MSE" = "#E74C3C", "MSE_BAY" = "#1F78B4")) +  # Set colors
  labs(
    x = "Number of Studies: n",
    y = "Mean Squared Error (MSE)",
    color = "Series"
  )


# Display the plot of MSE by sample size
x11()
print(p)
ggsave("Model_comp_FE_n50x500.pdf", plot = p, width = 8, height = 6, dpi = 300)

pooling(bg)
heterogeneity(bg)

### Partial Pooling -> 1.2.

for (n in seq(2, nstud, by = 2)) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s+3)
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0.1,max = 2*theta_true)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul[t])
    }
    # Frequentist
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    ## Using partial pooling from Bayesian bagger
    bg <- baggr(df_pooled, model="rubin", pooling = "partial")
    # Extract the Stan fit object
    fit_object <- bg$fit
    # Get the summary table from the Stan fit object
    fit_summary <- summary(fit_object)$summary
    # Convert it to a data frame for easier manipulation
    fit_summary_df <- as.data.frame(fit_summary)
    # Get row names containing "theta_k"
    theta_k_rows <- fit_summary_df[grep("theta_k", rownames(fit_summary)), ]
    theta_MSE <- (theta_k_rows$mean-1)^2
    BAY_MSE <- sum(unlist(theta_MSE))/n
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s)
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
  
  #bg_mse$n = (tau_columns[[n]]-1)^2
  #BAY_MSE[n] = sum(unlist(bg_mse))/n
}

# Plot the MSEs to compare
df <- data.frame(x=1:nstud, MLE_MSE, MSE_BAY)

df_long <- df %>%
  pivot_longer(cols = c(MLE_MSE, MSE_BAY), names_to = "Series", values_to = "Value")

# Plot with ggplot2
p2 <- ggplot(df_long, aes(x = x, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +  # Thicker lines for clarity
  geom_point(size = 2) +   # Distinct points for readability
  geom_text(aes(label = round(Value, 1)), 
            vjust = -0.7, hjust = 0.7, 
            size = 3.5, family = "Arial", color = "black",
            check_overlap = TRUE) +  # Prevent text overlap
  theme_minimal(base_family = "Arial") +  # Professional theme
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),  # Lighter gridlines
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("MLE_MSE" = "#E74C3C", "MSE_BAY" = "#1F78B4")) +  # Set colors
  labs(
    title = "MSE Comparison in a Fixed Effects Setting",
    x = "Number of Studies: n",
    y = "Mean Squared Error (MSE)",
    color = "Series"
  )


# Display the plot by sample size
x11()
print(p2)


### Full pooling -> 1.3.
for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s+3)
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0.1,max = 2*theta_true)
      theta_simul[t] = rnorm(1,mean = theta_true,sd = sigma_simul[t])
    }
    
    # Frequentist
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem#sum(theta_simul)
    
    # Bayesian
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    ## Using partial pooling from Bayesian bagger
    bg <- baggr(df_pooled, model="rubin", pooling = "full")
    # Extract the Stan fit object
    fit_object <- bg$fit
    # Get the summary table from the Stan fit object
    fit_summary <- summary(fit_object)$summary
    # Convert it to a data frame for easier manipulation
    fit_summary_df <- as.data.frame(fit_summary)
    # Get row names containing "theta_k"
    theta_k_rows <- fit_summary_df[grep("theta_k", rownames(fit_summary)), ]
    theta_MSE <- (theta_k_rows$mean-1)^2
    BAY_MSE <- sum(unlist(theta_MSE))/n
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s)
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
  
  #bg_mse$n = (tau_columns[[n]]-1)^2
  #BAY_MSE[n] = sum(unlist(bg_mse))/n
}

# Plot the MSEs to compare
df <- data.frame(x=1:nstud, MLE_MSE, MSE_BAY)

df_long <- df %>%
  pivot_longer(cols = c(MLE_MSE, MSE_BAY), names_to = "Series", values_to = "Value")

# Plot with ggplot2
p3 <- ggplot(df_long, aes(x = x, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +  # Thicker lines for clarity
  geom_point(size = 2) +   # Distinct points for readability
  geom_text(aes(label = round(Value, 1)), 
            vjust = -0.7, hjust = 0.7, 
            size = 3.5, family = "Arial", color = "black",
            check_overlap = TRUE) +  # Prevent text overlap
  theme_minimal(base_family = "Arial") +  # Professional theme
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),  # Lighter gridlines
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("MLE_MSE" = "#E74C3C", "MSE_BAY" = "#1F78B4")) +  # Set colors
  labs(
    title = "MSE Comparison in a Fixed Effects Setting",
    x = "Number of Studies: n",
    y = "Mean Squared Error (MSE)",
    color = "Series"
  )


# Display the plot
x11()
print(p3)


############ Graphs

bg2 <- baggr(df_pooled, model="rubin", pooling = "partial")
bg3 <- baggr(df_pooled, model="rubin", pooling = "none")
baggr_comparison <- baggr_compare("Partial pooling" = bg2, 
                                  "No pooling" = bg3, 
                                  plot = TRUE) 


PriorP1 <- baggr(df_pooled,ppd=TRUE)
PriorP2 <- baggr(schools, prior_hypermean = normal(300, 5), ppd=TRUE)
baggr_compare("Prior A, p.p.d."=PriorP1,
              "Prior B p.p.d."=PriorP2,
              compare = "effects", plot=TRUE)
#Isnt converging nicely -> High variance


help(baggr_compare)





############## Situation 4: DGM s.t. Regression #######

### Two groups

s=1         # Number of simulations
nstud = 20   # Maximum number of studies
theta_simul <- c()  #Simulate data set
theta_context <- c() #Simulate data set !NEW! Varying "true" ATE for each observation
sigma_simul <- c()  #Simulate data set
MLE_est <- c()    #Vector of simulated Frequentist weighted-pool Estimators
MLE_mean <- c()   #Monte Carlo average 
MLE_MSE <- c()    #Mean Squared Error of Frequentist approach
BAY_MSE <- c()    #Mean Squared Error of Bayesian approach
bg <- c()     #Bayesian Model
BAY_tau <- c()  
bg_mse <- c()   #Bayesian Mean Squared Error
MSE_hat <- 0
MSE_BAY <- c()
theta_true1 <- 1   #True Fixed Treatment Effect Group 1
theta_true2 <- 30  #True Fixed Treatment Effect Group 2
beta <- 3
covariates <-

for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s)
    for (t in 1:n) {    # Generator process (study-dependent)
      group_indicator <- sample(c(1, 2), 1)  # Randomly assign each t to a group 
      if (group_indicator == 1) {
        covariates[t] <- theta_true1
        dummy <- group indicator 
      } else {
        covariates[t] <- theta_true2
      }
      sigma_simul[t] = runif(1,min = 0.1,max = 0.5) #Variance
      theta_simul[t] = covariates[t]*beta1+dummy*beta2+rnorm(1,mean=0,sd=sigma_simul[t])
    }
    
    # Frequentist
    weighted_num = sum(theta_simul %*% sigma_simul)
    weighted_dem = sum(sigma_simul)
    MLE_est[i]<- weighted_num/weighted_dem #sum(theta_simul)
    
    # Bayesian
    df_pooled = data.frame("tau"= theta_simul,
                           "se" = sigma_simul)
    ## Using partial pooling from Bayesian bagger
    bg <- baggr(df_pooled, model="rubin", pooling="partial", iter=2000, chains=1)
    # Extract the Stan fit object
    fit_object <- bg$fit
    # Get the summary table from the Stan fit object
    fit_summary <- summary(fit_object)$summary
    # Convert it to a data frame for easier manipulation
    fit_summary_df <- as.data.frame(fit_summary)
    # Get row names containing "theta_k"
    theta_k_rows <- fit_summary_df[grep("theta_k", rownames(fit_summary)), ]
    theta_MSE <- (theta_k_rows$mean-theta_context)^2
    BAY_MSE <- sum(unlist(theta_MSE))/n
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s)
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
  
  #bg_mse$n = (tau_columns[[n]]-1)^2
  #BAY_MSE[n] = sum(unlist(bg_mse))/n
}

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
    title = "MSE Comparison in a Random Effects Setting",
    x = "Number of Studies (n)",
    y = "MSE"
  )

# Display the plot
x11()
print(p)

### Posterior distribution plot
effect_plot(bg)

### We can compare all three supported Bayesian Analysis in one graph
# my_baggr_comparison <- baggr_compare(df)
#plot(my_baggr_comparison) + 
#  ggtitle("8 schools: model comparison")




### Things to modify:
# Hyperpriors
# Data Generating Model: Variance assumptions etc. 
# Speed!
# Graphs for output
# For context: compare partial baggr with





