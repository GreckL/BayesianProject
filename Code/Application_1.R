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
total_consumption_per_cap_end_1_ETH = 0.239 # sd 0.068
Food_consumption_per_cap_end_1_ETH = 0.139 # 0.061
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_ETH = 0.347 # sd 0.074
Food_consumption_per_cap_end_2_ETH = 0.186 # sd 0.061

# Ghana P.67
# Treatment ITT - Endline 1
total_consumption_per_cap_end_1_GHN = 0.097 # sd 0.049
Food_consumption_per_cap_end_1_GHN = 0.065 # 0.044 (no star)
# Treatment ITT - Endline 2
total_consumption_per_cap_end_2_GHN = 0.136 # sd 0.050
Food_consumption_per_cap_end_2_GHN = 0.077 # sd 0.045 (one star)

# Honduras

