library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library("tidyverse")        # Data manipulation 
library("rstan")
library("compiler")
library("metafor")
rm(list=ls())

path <- "C:/Users/Grego/Documents/GitHub/BayesianProject"

################ Situation 2: DGM s.t. Random Effects Model with n ATEs ########

### Comparison: RE-Freq vs. Partial Pooling

s=50       # Number of simulations
nstud = 40   # Maximum number of studies
MLE_est <- c()    #Vector of simulated Frequentist weighted-pool Estimators
BAY_MSE <- c()    #Mean Squared Error of Bayesian approach
bg <- c()     #Bayesian Model
BAY_tau <- c()  
theta_true <- 1   #True Average Treatment Effect
RMA <- numeric(nstud)
MSE_hat <- numeric(nstud)

custom_prior  <- list(hypermean = normal(0,8^2), hypersd = uniform(0,10*0.25))

set.seed(123)            # For reproducibility (if desired)

###############################################################################
# Main Loop: Vary the Number of Studies (n), Then Do Multiple Simulations (s)
###############################################################################
for (n in seq(2, nstud, by = 1)) {
  
  for (i in 1:s) {  

    sigma_simul    <- numeric(n)  # Standard deviations
    theta_context  <- numeric(n)  # "True" effects for each study
    theta_simul    <- numeric(n)  # Observed effects for each study
    
    for (t in 1:n) {
      sigma_simul[t]    <- rhalfnorm(1, theta = 5)
      theta_context[t]  <- rnorm(1, mean = theta_true, sd = 1)
      theta_simul[t]    <- rnorm(1, mean = theta_context[t], 
                                 sd   = sigma_simul[t])
    }
   
    # FREQUENTIST META-ANALYSIS (metafor)
    res       <- rma.uni(yi = theta_simul, sei = sigma_simul, method = "REML")
    blup_df   <- blup(res)
    blup_preds <- blup_df$pred
    
    # MSE for the frequentist BLUPs relative to the true effects
    sq_errors_freq <- (blup_preds - theta_context)^2
    MSE_freq       <- mean(sq_errors_freq)
    
    # Incrementally update the Monte Carlo average for frequentist approach
    RMA[n] <- RMA[n] * ((i - 1) / i) + MSE_freq * (1 / i)
    
    # BAYESIAN META-ANALYSIS (baggr)
    df_pooled <- data.frame(
      tau = theta_simul,    # Observed effects
      se  = sigma_simul     # Their standard errors
    )
    bg <- baggr(df_pooled, model = "rubin", pooling = "partial", iter = 2000,
                prior = custom_prior)
    # Extract the Stan fit & summary
    fit_object   <- bg$fit
    fit_summary  <- summary(fit_object)$summary
    fit_df       <- as.data.frame(fit_summary)
    
    # Identify rows corresponding to the study-specific parameters in baggr
    idx_thetas   <- grep("theta_k", rownames(fit_df))
    theta_k_rows <- fit_df[idx_thetas, , drop = FALSE]
    
    # Check dimension: we expect n rows in 'theta_k_rows'
    if (nrow(theta_k_rows) != n) {
      warning(sprintf("Number of 'theta_k' rows (%d) != n (%d)!", 
                      nrow(theta_k_rows), n))
    }
    # MSE for the Bayesian partial pooling
    theta_MSE <- (theta_k_rows$mean - theta_context)^2
    BAY_MSE   <- mean(theta_MSE)  # or sum(...) / n, same result

    # Incrementally update the Monte Carlo average for Bayesian approach
    MSE_hat[n] <- MSE_hat[n] * ((i - 1) / i) + BAY_MSE * (1 / i)
  }
  # Print progress
  print(sprintf("Finished n = %d", n))
}

RMA        # Frequentist MSE for each n
MSE_hat     # Bayesian MSE for each n

df_plotImportant <- data.frame(
  n = 1:nstud,
  RMA = RMA,
  BAYES = MSE_hat
)


###############################################################################
# 2) Reshape (pivot) your data so that the value column is "RMA"
###############################################################################
df_plotImportant_filtered <- df_plotImportant %>%
  filter(n != 1)

# Reshape (pivot) the data
dflong_Imp <- df_plotImportant_filtered %>%
  pivot_longer(
    cols = c(RMA, BAYES),      # Corrected column names to match your DataFrame
    names_to = "Series",       # Column to store series names
    values_to = "Value"        # Column to store corresponding values
  )

# Create the ggplot
p <- ggplot(dflong_Imp, aes(x = n, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +          # Thicker lines
  geom_point(size = 2) +           # Points
  theme_minimal(base_family = "sans") +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c("RMA" = "#E74C3C", 
                                "BAYES" = "#1F78B4")) +
  labs(
    x = "Number of Studies: n",
    y = "Mean Squared Error (MSE)",
    color = "Series"
  )

# Display the plot
print(p)
ggsave("Model_comp_RE_n50x500.pdf", plot = p, width = 8, height = 6, dpi = 300)



########### Priors
custom_prior  <- list(hypermean = normal(0,8^2), hypersd = uniform(0,10*0.25))

PriorP1 <- baggr(df_pooled,prior=custom_prior,ppd=TRUE)
PriorP2 <- baggr(df_pooled, prior_hypermean = normal(1, 20), 
                 prior_hypersd = uniform(0,10*0.25), ppd=TRUE)
PriorComp=baggr_compare("Prior A"=PriorP1,
              "Prior B"=PriorP2,
              compare = "effects", plot=TRUE)






















############## Situation 3: DGM s.t. Random Effects Model with group ATE #######

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

for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s)
    for (t in 1:n) {    # Generator process (study-dependent)
      group_indicator <- sample(c(1, 2), 1)  # Randomly assign each t to a group 
      if (group_indicator == 1) {
        theta_context[t] <- theta_true1
      } else {
        theta_context[t] <- theta_true2
      }
      sigma_simul[t] = runif(1,min = 0.1,max = 0.5) #Variance
      theta_simul[t] = rnorm(1,mean = theta_context[t],sd = sigma_simul[t]) #ATE
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