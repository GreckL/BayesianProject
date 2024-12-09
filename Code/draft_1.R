###
install.packages("baggr")
library("bayesmeta")
library("baggr")

### 1 - Simulating data ###
n = 20
tau_simul = runif(n,min = -1,max = 1)
sigma_simul = runif(n,min = .5,max = 5)



df_pooled <- data.frame("tau" = tau_simul,
                        "se"  = sigma_simul)
bg <- baggr(df_pooled, pooling = "full")