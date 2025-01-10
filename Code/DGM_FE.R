### Bayesian Meta Analysis
#install.packages("baggr")
#install.packages("bayesmeta")
#install.packages("bayesplot")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library("tidyverse")        # Data manipulation 
library("rstan")


rm(list=ls())
set.seed(12345)

########## Situation 1: DGM s.t. Fixed Effects and no Covariates ###############

# Parameters
s=50         # Number of simulations
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
df <- df %>%
  mutate(MSE_MLE = MLE_MSE)
dflong <- df %>%
  pivot_longer(cols = c(MSE_MLE, MSE_BAY), names_to = "Series", values_to = "Value")



# Plot with ggplot2
p <- ggplot(dflong, aes(x = x, y = Value, color = Series, group = Series)) +
  geom_line(size = 1.2) +   # Thicker lines
  geom_point(size = 2) +    # Points
  # Remove or comment out the geom_text layer:
  # geom_text(aes(label = round(Value, 1)), 
  #           vjust = -0.7, hjust = 0.7, 
  #           size = 3.5, family = "sans", color = "black",
  #           check_overlap = TRUE) +
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
  scale_color_manual(values = c("MSE_MLE" = "#E74C3C", "MSE_BAY" = "#1F78B4")) +
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
