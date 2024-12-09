### Load libraries
library("bayesmeta")
library("baggr")

### 1 - Simulating data ###
n = 20
tau_simul = runif(n,min = -1,max = 1)
sigma_simul = runif(n,min = -.5,max = .5)




library(baggr)
df_pooled <- data.frame("tau" = c(28,8,-3,7,-1,1,18,12),
                        "se"  = c(15,10,16,11,9,11,10,18))
bg <- baggr(df_pooled, pooling = "partial")