rm(list=ls())

source("CNP_lake_functions.R")


# Step 0: Replicating Tyrell results ----
p=Get_classical_param_lake()
ini=Get_initial_values(p)
ini[1:3]=0
p[c("ID","IP","IN")]=3
d=Compute_ode(ini,p,julia = F,n_time = 100,solver = "ode45")


