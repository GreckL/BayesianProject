### Bayesian Meta Analysis - Application

library("bayesmeta")
library("baggr")
library("bayesplot")
library("ggplot2")        # For fancy plotting 
library("tidyverse")        # Data manipulation 
library("rstan")
library("compiler")

rm(list=ls())
set.seed(12345)

### 1 - Simulating data ###


################## Situation 1: DGM s.t. Fixed Effects ########################
## Total meta analysis result
# Table 3: Per Capita Consumption P. 57
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1 = 4.55
Food_consumption_per_cap_end_1 = 3.87
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2 = 3.36
Food_consumption_per_cap_end_2 = 2.62

# Each country result
# Ethiopia P.66
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_ETH = 0.239 # sd 0.068
Food_consumption_per_cap_end_1_std_ETH = 0.139 # 0.061
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_ETH = 0.347 # sd 0.074
Food_consumption_per_cap_end_2_std_ETH = 0.186 # sd 0.061

# Ghana P.67
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_GHN = 0.097 # sd 0.049
Food_consumption_per_cap_end_1_std_GHN = 0.065 # sd 0.044 (no star)
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_GHN = 0.136 # sd 0.050
Food_consumption_per_cap_end_2_std_GHN = 0.077 # sd 0.045 (one star)

# Honduras P.68
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_HON = 0.011 # sd 0.049 (no star)
Food_consumption_per_cap_end_1_std_HON = 0.136 # sd 0.047 
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_HON = -0.070 # sd 0.049 (no star)
Food_consumption_per_cap_end_2_std_HON = 0.088 # sd 0.051 (one star)

# India P.69
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_IND = 0.296 # sd 0.080 
Food_consumption_per_cap_end_1_std_IND = 0.238 # sd 0.068 
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_IND = 0.228 # sd 0.066
Food_consumption_per_cap_end_2_std_IND = 0.278 # sd 0.063

# Pakistan P.70
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_PAK = 0.171 # sd 0.064 
Food_consumption_per_cap_end_1_std_PAK = 0.117 # sd 0.056 (two stars) 
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_PAK = 0.117 # sd 0.067 (one star)
Food_consumption_per_cap_end_2_std_PAK = 0.058 # sd 0.059 (no star)

# Peru P.71
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_std_PER = 0.048 # sd 0.048 (no star)
Food_consumption_per_cap_end_1_std_PER = 0.020 # sd 0.052 (no star)
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_std_PER = 0.096 # sd 0.056 (one star)
Food_consumption_per_cap_end_2_std_PER = 0.064 # sd 0.045 (no star)

## Per Capita Consumption Per Country P.72-73
# Ethiopia 
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_ETH = 6.83 # sd 1.93
sd_total_consumption_per_cap_end_1_ETH = 1.93
Food_consumption_per_cap_end_1_ETH = 2.42 # sd 1.48 (no stars)
sd_Food_consumption_per_cap_end_1_ETH = 1.48
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_ETH = 7.37 # sd 1.58
sd_total_consumption_per_cap_end_2_ETH = 1.58
Food_consumption_per_cap_end_2_ETH = 2.02 # sd 0.86
sd_Food_consumption_per_cap_end_2_ETH = 0.86

# Ghana 
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_GHN = 2.82 # sd 1.42 (two stars)
sd_total_consumption_per_cap_end_1_GHN = 1.42
Food_consumption_per_cap_end_1_GHN = 2.18 # sd 1.04 (two stars)
sd_Food_consumption_per_cap_end_1_GHN = 1.04
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_GHN = 3.22 # sd 1.19
sd_total_consumption_per_cap_end_2_GHN = 1.19
Food_consumption_per_cap_end_2_GHN = 2.41 # sd 0.81
sd_Food_consumption_per_cap_end_2_GHN = 0.81

# Honduras P.68
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_HON = 0.61 # sd 2.70 (no star)
sd_total_consumption_per_cap_end_1_HON = 2.70
Food_consumption_per_cap_end_1_HON = 3.45 # sd 1.51 (two stars) 
sd_Food_consumption_per_cap_end_1_HON = 1.51
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_HON = -4.45 # sd 3.11 (no star)
sd_total_consumption_per_cap_end_2_HON = 3.11
Food_consumption_per_cap_end_2_HON = -0.60 # sd 1.28 (no star)
sd_Food_consumption_per_cap_end_2_HON = 1.28

# India P.69
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_IND = 6.51 # sd 1.75 
sd_total_consumption_per_cap_end_1_IND = 1.75
Food_consumption_per_cap_end_1_IND = 4.96 # sd 1.17
sd_Food_consumption_per_cap_end_1_IND = 1.17
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_IND = 6.18 # sd 1.78
sd_total_consumption_per_cap_end_2_IND = 1.78
Food_consumption_per_cap_end_2_IND = 3.02 # sd 1.11
sd_Food_consumption_per_cap_end_2_IND = 1.11

# Pakistan P.70
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_PAK = 8.86 # sd 3.31 
sd_total_consumption_per_cap_end_1_PAK = 3.31
Food_consumption_per_cap_end_1_PAK = 3.92 # sd 1.78 (two stars) 
sd_Food_consumption_per_cap_end_1_PAK = 1.78
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_PAK = 5.98 # sd 3.40 (one star)
sd_total_consumption_per_cap_end_2_PAK = 3.40
Food_consumption_per_cap_end_2_PAK = 2.91 # sd 1.80 (no star)
sd_Food_consumption_per_cap_end_2_PAK = 1.80

# Peru P.71
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_PER = 4.18 # sd 4.14 (no star)
sd_total_consumption_per_cap_end_1_PER = 4.14
Food_consumption_per_cap_end_1_PER = 5.70 # sd 3.22 (one star)
sd_Food_consumption_per_cap_end_1_PER = 3.22
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_PER = 6.18 # sd 3.59 (one star)
sd_total_consumption_per_cap_end_2_PER = 3.59
Food_consumption_per_cap_end_2_PER = 6.25 # sd 2.58 (two star)
sd_Food_consumption_per_cap_end_2_PER = 2.58

#build data frame
#taus
total_cons_endline1 = c(total_consumption_per_cap_end_1_ETH, total_consumption_per_cap_end_1_GHN,
                        total_consumption_per_cap_end_1_HON, total_consumption_per_cap_end_1_IND,
                        total_consumption_per_cap_end_1_PAK, total_consumption_per_cap_end_1_PER)
total_food_endline1 = c(Food_consumption_per_cap_end_1_ETH, Food_consumption_per_cap_end_1_GHN,
                        Food_consumption_per_cap_end_1_HON, Food_consumption_per_cap_end_1_IND,
                        Food_consumption_per_cap_end_1_PAK, Food_consumption_per_cap_end_1_PER)
total_cons_endline2 = c(total_consumption_per_cap_end_2_ETH, total_consumption_per_cap_end_2_GHN,
                        total_consumption_per_cap_end_2_HON, total_consumption_per_cap_end_2_IND,
                        total_consumption_per_cap_end_2_PAK, total_consumption_per_cap_end_2_PER)
total_food_endline2 = c(Food_consumption_per_cap_end_2_ETH, Food_consumption_per_cap_end_2_GHN,
                        Food_consumption_per_cap_end_2_HON, Food_consumption_per_cap_end_2_IND,
                        Food_consumption_per_cap_end_2_PAK, Food_consumption_per_cap_end_2_PER)
#SDs
sd_total_cons_endline1 = c(sd_total_consumption_per_cap_end_1_ETH, sd_total_consumption_per_cap_end_1_GHN,
                           sd_total_consumption_per_cap_end_1_HON, sd_total_consumption_per_cap_end_1_IND,
                           sd_total_consumption_per_cap_end_1_PAK, sd_total_consumption_per_cap_end_1_PER)
sd_total_food_endline1 = c(sd_Food_consumption_per_cap_end_1_ETH, sd_Food_consumption_per_cap_end_1_GHN,
                        sd_Food_consumption_per_cap_end_1_HON, sd_Food_consumption_per_cap_end_1_IND,
                        sd_Food_consumption_per_cap_end_1_PAK, sd_Food_consumption_per_cap_end_1_PER)
sd_total_cons_endline2 = c(sd_total_consumption_per_cap_end_2_ETH, sd_total_consumption_per_cap_end_2_GHN,
                        sd_total_consumption_per_cap_end_2_HON, sd_total_consumption_per_cap_end_2_IND,
                        sd_total_consumption_per_cap_end_2_PAK, sd_total_consumption_per_cap_end_2_PER)
sd_total_food_endline2 = c(sd_Food_consumption_per_cap_end_2_ETH, sd_Food_consumption_per_cap_end_2_GHN,
                        sd_Food_consumption_per_cap_end_2_HON, sd_Food_consumption_per_cap_end_2_IND,
                        sd_Food_consumption_per_cap_end_2_PAK, sd_Food_consumption_per_cap_end_2_PER)

# Bayesian approach
#cons endline 1
df_pooled_cons_endline1 = data.frame("tau"= total_cons_endline1,
                       "se" = sd_total_cons_endline1)
mean_bay_cons_endline1 = baggr(df_pooled_cons_endline1, model="rubin", pooling = "partial")
bay_outcome_cons_endline1 = pooling(mean_bay_cons_endline1)
#food endline 1
df_pooled_food_endline1 = data.frame("tau"= total_food_endline1,
                                     "se" = sd_total_food_endline1)
mean_bay_food_endline1 = baggr(df_pooled_food_endline1, model="rubin", pooling = "partial")
bay_outcome_food_endline1 = pooling(mean_bay_food_endline1)
#cons endline 2
df_pooled_cons_endline2 = data.frame("tau"= total_cons_endline2,
                                     "se" = sd_total_cons_endline2)
mean_bay_cons_endline2 = baggr(df_pooled_cons_endline2, model="rubin", pooling = "partial")
bay_outcome_cons_endline2 = pooling(mean_bay_cons_endline2)
#food endline 2
df_pooled_food_endline2 = data.frame("tau"= total_food_endline2,
                                     "se" = sd_total_food_endline2)
mean_bay_food_endline2 = baggr(df_pooled_food_endline2, model="rubin", pooling = "partial")
bay_outcome_food_endline2 = pooling(mean_bay_food_endline2)



# Extract the Stan fit object
fit_object <- mean_bay_cons_endline1$fit
# Get the summary table from the Stan fit object
fit_summary <- summary(fit_object)$summary
# Convert it to a data frame for easier manipulation
fit_summary_df <- as.data.frame(fit_summary)
################ Situation 2: DGM s.t. Random Effects Model with n ATEs ########

### Comparison: FE-Freq vs. Partial Pooling

s=1         # Number of simulations
nstud = 6   # Number of studies
theta_simul <- c()  #Simulate data set (in this case we would need the microdata)
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
theta_true <- 1   #True Average Treatment Effect

for (n in 2:nstud) {    # Number of studies (meta analysis so start from 2 up to nstud)
  for (i in 1:s) {      # Number of simulations
    set.seed(s)
    for (t in 1:n) {    # Generator process (study-dependent)
      sigma_simul[t] = runif(1,min = 0,max = 0.5) #Variance
      theta_context[t] = rnorm(1, mean=theta_true, sd = 1) 
      theta_simul[t] = rnorm(1,mean = theta_wacontext[t],sd = sigma_simul[t]) #ATE
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
    
    theta_MSE <- (theta_k_rows$mean-theta_context)^2
    BAY_MSE <- sum(unlist(theta_MSE))/n
    MSE_hat = MSE_hat*((s-1)/s) + BAY_MSE*(1/s) #Average over simulations
    #*BAY_Nstud_MSE[i] <- sum(BAY_MSE)/n
    # Print the theta_k rows
    #print(theta_k_rows)
  }
  MLE_mean[n] <- sum(MLE_est)/s #Average MLE
  MLE_MSE[n] <- (MLE_mean[n]-theta_true)^2 #MSE of average MLE -> But we want Average MSE
  
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