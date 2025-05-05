rm(list=ls())
source("CNP_lake_functions.R")

Run_with_julia=T



# ----------------------------- Step 1: Fixers and decomposers ----
## >> 1) Gradient IP, ID, IN, fixed stoichio ----

d2=d_indirect=tibble()
n_=100
param_list=rbind(expand.grid(IN=seq(1,50,length.out=n_),
                             IP=c(5),
                             ID=c(5))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=5,
                             IP=seq(1,50,length.out=n_),
                             ID=c(5))%>%
                   add_column(., Simulation_ID=2),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_))%>%
                   add_column(., Simulation_ID=3))

p=Get_classical_param_lake()
ini=Get_initial_values(p)
ini[c(7:9)]=0

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  p$same_stoichio=T
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  if (max(d$Time)!=10000){print("Not_converged")}
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,]))
  
  print(x)
}

write.table(d2,"./data/Simulations/Simulation_fixers_and_decomposers_fixed.csv",sep=";")

## >> 2) Gradient ID quality ----

d2=d_indirect=tibble()
n_=40
param_list=rbind(expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_),
                             beta_A=c(.002,.012),
                             alpha_A=seq(.02,.14,length.out=2))%>%
                   add_column(., Simulation_ID="Alpha_vary"),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_),
                             beta_A=seq(.002,.012,length.out=2),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID="Beta_vary"))

p=Get_classical_param_lake()
ini=Get_initial_values(p)
ini[c(7:9)]=0

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  p$same_stoichio=F
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$alpha_allo=param_list$alpha_A[x]
  p$beta_allo=param_list$beta_A[x]
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  if (max(d$Time)!=10000){print("Not_converged")}
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,]))
  
  print(x)
}

write.table(d2,"./data/Simulations/Simulation_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")

# ----------------------------- Step 2: Non-Fixers and decomposers ----
## >> 1) Gradient IP, ID, IN, fixed stoichio ----

d2=d_indirect=tibble()
n_=100
param_list=rbind(expand.grid(IN=seq(1,50,length.out=n_),
                             IP=c(5),
                             ID=c(5))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=5,
                             IP=seq(1,50,length.out=n_),
                             ID=c(5))%>%
                   add_column(., Simulation_ID=2),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_))%>%
                   add_column(., Simulation_ID=3))

p=Get_classical_param_lake()
ini=Get_initial_values(p)
ini[c(4:6)]=0

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  p$same_stoichio=T
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  if (max(d$Time)!=10000){print("Not_converged")}
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,]))
  
  print(x)
}

write.table(d2,"./data/Simulations/Simulation_Non_fixers_and_decomposers_fixed.csv",sep=";")

## >> 2) Gradient ID quality ----

d2=d_indirect=tibble()
n_=40
param_list=rbind(expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_),
                             beta_A=c(.002,.012),
                             alpha_A=seq(.02,.14,length.out=2))%>%
                   add_column(., Simulation_ID="Alpha_vary"),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_),
                             beta_A=seq(.002,.012,length.out=2),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID="Beta_vary"))

p=Get_classical_param_lake()
ini=Get_initial_values(p)
ini[c(4:6)]=0

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  p$same_stoichio=F
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$alpha_allo=param_list$alpha_A[x]
  p$beta_allo=param_list$beta_A[x]
  
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  if (max(d$Time)!=10000){print("Not_converged")}
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,]))
  
  print(x)
}

write.table(d2,"./data/Simulations/Simulation_Non_fixers_and_decomposers_2D_quality_quantity.csv",sep=";")



# ----------------------------- Step 3: CNP allochtonous flows ----
## >> 1/ Increasing IN, IP, ID ----

d2=d_alone=tibble()
n_=100
compute_alone=T
param_list=rbind(expand.grid(IN=seq(1,10,length.out=n_),
                             IP=c(5),
                             ID=c(25),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=5,
                             IP=seq(.5,5,length.out=n_),
                             ID=c(5),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=2),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,100,length.out=n_),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=3))

p=Get_classical_param_lake()

for (x in 1:nrow(param_list)){  
  ini=Get_initial_values(p)
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$beta_allo=param_list$beta_A[x]
  p$alpha_allo=param_list$alpha_A[x]
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)!=4000,F,T)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Get_CNP(d,p),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_,
                    p))
  
  if (compute_alone){
    save_ini=ini
    
    ini[4:9]=0
    
    #Only decompo
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="Decomposers"))
    
    #Only fixers
    
    ini=save_ini
    ini[c(1:3,7:9)]=0
    
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="Fixers"))
    
    #Only non fixers
    
    ini=save_ini
    ini[c(1:6)]=0
    
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="Non_fixers"))
    
    
    # F and decompo
    ini=save_ini
    ini[c(7:9)]=0
    
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="F_and_B"))
    
    
    # NF and decompo
    ini=save_ini
    ini[c(4:6)]=0
    
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="NF_and_B"))
    
    # F and NF
    ini=save_ini
    ini[c(1:3)]=0
    
    d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    converged_=ifelse(max(d$Time)!=4000,F,T)
    
    Eq=Get_equilibrium(d,p)
    d_alone=rbind(d_alone,cbind(Get_equilibrium(d,p),
                                Get_limitation(d,p),
                                Get_CNP(d,p),
                                Inflow="IN",
                                Get_phi_functions(d,p),
                                param_list[x,],
                                converged_,
                                Which_species_alone="NF_and_F"))
    
  }
  
  
  print(x)
}

write.table(d2,"./data/Simulation_allochtonous.csv",sep=";")
write.table(d_alone,"./data/Simulation_allochtonous_species_alone.csv",sep=";")


## >> 2/ Increasing alpha_A, beta_A ----

d2=d_alone=tibble()
n_=40
compute_alone=T
param_list=rbind(expand.grid(IN=c(5),
                             IP=c(5),
                             ID=seq(1,65,length.out=n_),
                             beta_A=seq(.002,.012,length.out=n_),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=c(5),
                             IP=c(5),
                             ID=seq(1,65,length.out=n_),
                             beta_A=c(.002,.012),
                             alpha_A=seq(.02,.14,length.out=n_))%>%
                   add_column(., Simulation_ID=2))

p=Get_classical_param_lake()

for (x in 1:nrow(param_list)){  
  ini=Get_initial_values(p)
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$beta_allo=param_list$beta_A[x]
  p$alpha_allo=param_list$alpha_A[x]
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)!=4000,F,T)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="A_to_I"),
                    Get_CNP(d,p),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_,
                    p))
  


write.table(d2,"./data/Simulation_allochtonous_alpha_beta_A.csv",sep=";")



# ----------------------------- Step 4: Other indirect effects ----

d2=d_alone=tibble()
n_=100
compute_alone=T
param_list=rbind(expand.grid(IN=seq(1,10,length.out=n_),
                             IP=c(5),
                             ID=c(25),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=5,
                             IP=seq(.05,10,length.out=n_),
                             ID=c(5),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=2),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,100,length.out=n_),
                             beta_A=c(.002,.012),
                             alpha_A=c(.02,.14))%>%
                   add_column(., Simulation_ID=3))

p=Get_classical_param_lake()

for (x in 1:nrow(param_list)){  
  ini=Get_initial_values(p)
  
  p$functional_response_phyto=1
  
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$beta_allo=param_list$beta_A[x]
  p$alpha_allo=param_list$alpha_A[x]
  
  d=Compute_ode(ini,p,julia = Run_with_julia,n_time = 4000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)!=4000,F,T)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Indirect_effects(d,p)%>%filter(., Type=="I_to_A"),
                    Get_CNP(d,p),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_,
                    p))
  
  print(x)
}

write.table(d2,"./data/Simulation_allochtonous_I_to_A.csv",sep=";")

