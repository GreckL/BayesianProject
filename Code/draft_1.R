### Load libraries
library("bayesmeta")
library("baggr")

### 1 - Simulating data ###
n = 20
tau_simul = runif(n,min = -1,max = 1)
sigma_simul = runif(n,min = -.5,max = .5)
