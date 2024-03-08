rm(list=ls())

source("CNP_lake_functions.R")
# for (k in 1:length(p)){assign(names(p)[k],p[[k]])}
# for (k in 1:length(ini)){assign(names(ini)[k],ini[k])}



p=Get_classical_param_lake()
ini=Get_initial_values(p)

p$ID=50
p$IN=10
p$IP=3
p$lP=2
p$lD=1
p$dB=.2
p[c("sB","sO","sF")]=0
p$DC_=F
d=Compute_ode(ini,p,julia = F,n_time = 2000)
Plot_D(Add_CNP(d,p))
Plot_net_effects(d,p)
Get_limitation(d,p)
Get_equilibrium(d,p)
Get_CNP(d,p)
Get_phi_functions(d,p)


# Step 1: CNP allochtonous flows ----
## >> 1/ Increasing IN, IP, ID ----

d2=d_indirect=tibble()
n_=100
param_list=rbind(expand.grid(IN=seq(1,50,length.out=n_),
                             IP=c(5),
                             ID=c(25))%>%
                   add_column(., Simulation_ID=1),
                 expand.grid(IN=5,
                             IP=seq(.05,10,length.out=n_),
                             ID=c(5))%>%
                   add_column(., Simulation_ID=2),
                 expand.grid(IN=5,
                             IP=c(5),
                             ID=seq(1,50,length.out=n_))%>%
                   add_column(., Simulation_ID=3))

p=Get_classical_param_lake()
ini=Get_initial_values(p)

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  
  p[c("sB","sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  
  
  d=Compute_ode(ini,p,julia = T,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  if (max(d$Time)!=10000){print("Not_converged")}
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Inflow="IN",
                    Get_phi_functions(d,p),
                    param_list[x,]))

  d_indirect=rbind(d_indirect,
                   cbind(Indirect_effects(d,p),
                         Get_limitation(d,p),
                         Get_CNP(d,p),
                         Inflow="IN",
                         param_list[x,]
                         ))
  print(x)
}

write.table(d2,"../Table/Simulation_allochtonous.csv",sep=";")
write.table(d_indirect,"../Table/Simulation_allochtonous_indirect.csv",sep=";")



## >> 2/ 2D increase of parameters ----

d2=d_indirect=tibble()
n_=50
param_list=expand.grid(IN=seq(1,50,length.out=n_),
                       ID=seq(1,50,length.out=n_),
                       IP=c(1.2,5,10))

p=Get_classical_param_lake()
ini=Get_initial_values(p)


for (x in 1:nrow(param_list)){
  
  p=Get_classical_param_lake()
  ini=Get_initial_values(p)
  
  p$IN=param_list$IN[x]
  p$IP=param_list$IP[x]
  p$ID=param_list$ID[x]
  p$functional_response_phyto=1
  p[c("sB","sF","sO")]=0.01
  
  d=Compute_ode(ini,p,julia = T,n_time = 5000)
  converged_=ifelse(max(d$Time)==5000,T,F)
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    param_list[x,],
                    converged_))
  

  d_indirect=rbind(d_indirect,
                   cbind(Indirect_effects(d,p),
                         Get_limitation(d,p),
                         Get_CNP(d,p),
                         param_list[x,],
                         converged_))
  print(x)
  
}

write.table(d2,"../Table/Simulation_allochtonous_2D.csv",sep=";")
write.table(d_indirect,"../Table/Simulation_allochtonous_indirect_2D.csv",sep=";")





# Step 2: CNP allochtonous flows ratio ----

d2=d_indirect=tibble()
n_=100
param_list=rbind(
  expand.grid(alpha_allo=.05,
              beta_allo=seq(.001,.015,length.out=n_),
              ID=c(1,10,50),
              Ratio_changed="PC"),
  expand.grid(alpha_allo=seq(.01,.15,length.out=n_),
              beta_allo=.005,
              ID=c(1,10,50),
              Ratio_changed="NC")
)

for (x in 1:nrow(param_list)){
  
  p=Get_classical_param_lake()
  
  p$functional_response_phyto=1
  p[c("sB","sF","sO")]=0.1
  
  p$ID=param_list$ID[x]
  p$functional_response_phyto=1
  p$beta_allo=param_list$beta_allo[x]
  p$alpha_allo=param_list$alpha_allo[x]
  
  d=Compute_ode(ini,p,julia = T,n_time = 10000)
  
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    param_list[x,]))
  
  
  d_indirect=rbind(d_indirect,
                   cbind(Indirect_effects(d,p),
                         Get_limitation(d,p),
                         Get_CNP(d,p),
                         param_list[x,]
                   ))
  print(x)
  
}

write.table(d2,"../Table/Simulation_allochtonous_ratios.csv",sep=";")
write.table(d_indirect,"../Table/Simulation_allochtonous_indirect_ratios.csv",sep=";")



# Step 3: Increasing self-regulation all ----

d2=tibble()
n_=100
param_list=expand.grid(IN=c(5,20),
                       IP=c(5),
                       ID=c(5,40),
                       self=seq(.01,1,length.out=100))

p=Get_classical_param_lake()
ini=Get_initial_values(p)

for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  
  p[c("sB","sF","sO")]=param_list$self[x]
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  
  
  d=Compute_ode(ini,p,julia = T,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)==10000,T,F)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Inflow="s",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_))
  
  print(x)
}

write.table(d2,"../Table/Simulation_allochtonous_self_regulation.csv",sep=";")

# Step 4: Increasing self-regulation, each ----

d2=tibble()
n_=100
param_list=expand.grid(IN=c(5,20),
                       IP=c(5),
                       ID=c(5,40),
                       self=seq(.01,1,length.out=100),
                       organism=c("B","F","O"))


for (x in 1:nrow(param_list)){  
  
  # if (type_flows=="IP") p$IN=50
  p=Get_classical_param_lake()
  ini=Get_initial_values(p)
  
  p$functional_response_phyto=1
  
  if (param_list$organism[x]=="B"){
    p[c("sB")]=param_list$self[x]
  }else if (param_list$organism[x]=="F"){
    p[c("sF")]=param_list$self[x]
  }else{
    p[c("sO")]=param_list$self[x]
  }
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  
  
  d=Compute_ode(ini,p,julia = T,n_time = 10000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)==10000,T,F)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Inflow="s",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_))
  
  print(x)
}

write.table(d2,"../Table/Simulation_allochtonous_self_regulation_each.csv",sep=";")



# Step 5: Stoichiometric ratio decomposers & self-regulation decomposers ----

d2=tibble()
param_list=expand.grid(IN=c(5,20),
                       IP=c(5),
                       ID=c(5,40),
                       alphaB=seq(.15,.25,length.out=30),
                       sB=seq(.01,1,length.out=30))


for (x in 1:nrow(param_list)){  
  p=Get_classical_param_lake()
  ini=Get_initial_values(p)
  
  # if (type_flows=="IP") p$IN=50
  
  p$functional_response_phyto=1
  
  p[c("sF","sO")]=0.1
  p$IN=param_list$IN[x]
  p$ID=param_list$ID[x]
  p$IP=param_list$IP[x]
  p$alpha_B=param_list$alphaB[x]
  p$sB=param_list$sB[x]
  
  d=Compute_ode(ini,p,julia = T,n_time = 5000)
  #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
  converged_=ifelse(max(d$Time)==5000,T,F)
  
  Eq=Get_equilibrium(d,p)
  d2=rbind(d2,cbind(Get_equilibrium(d,p),
                    Get_limitation(d,p),
                    Get_CNP(d,p),
                    Inflow="s",
                    Get_phi_functions(d,p),
                    param_list[x,],
                    converged_))
  
  print(x)
}

write.table(d2,"../Table/Simulation_traits_decomposers.csv",sep=";")

# Step 6: Random parameters ----


Run_random_param_function=function(ID){

  source("CNP_lake_functions.R")

  Random_param_list=matrix(0,nrow = 8,ncol = 2) 
  rownames(Random_param_list)=c('ID',"IP", "IN", "sO","sF","sB","alpha_B","beta_B")
  Random_param_list[,1]=c(3,3,3,.01,.01,.01,.12,.015)
  Random_param_list[,2]=c(50,50,50,.1,.1,.1,.25,.025)
  Random_param_list=as.data.frame(Latinhyper(Random_param_list, 500))
  
  d2=tibble()
  
  for (x in 1:nrow(Random_param_list)){  
    
    p=Get_classical_param_lake()
    ini=Get_initial_values(p)
    
    # if (type_flows=="IP") p$IN=50
    
    p$functional_response_phyto=1
    
    p$IN=Random_param_list$IN[x]
    p$ID=Random_param_list$ID[x]
    p$IP=Random_param_list$IP[x]
  
    p$sB=Random_param_list$sB[x]
    p$sO=Random_param_list$sO[x]
    p$sF=Random_param_list$sF[x]
  
    p$alpha_B=Random_param_list$alpha_B[x]
    p$beta_B=Random_param_list$beta_B[x]
    
    
    
    d=Compute_ode(ini,p,julia = T,n_time = 5000)
    #d=Compute_ode(ini,p,julia = F,n_time = 5000,solver = "lsoda")
    if (max(d$Time)!=5000){print("Not_converged")}
    
    Eq=Get_equilibrium(d,p)
    d2=rbind(d2,cbind(Get_equilibrium(d,p),
                      Get_limitation(d,p),
                      Get_CNP(d,p),
                      Get_phi_functions(d,p),
                      Random_param_list[x,]))
    
    print(x)
  }

  write.table(d2,paste0("../Table/Random/Simulation_allochtonous_random_param_",ID,".csv"),sep=";")

}

library(parallel)
mclapply(1:60,Run_random_param_function,mc.cores = 60)

#Aggregating data

d=tibble()
list_f=list.files("../Table/Random",pattern = ".csv")
for (k in list_f){
  d=rbind(d,read.table(paste0("../Table/Random/",k),sep=";"))
}

d$Deviation_Redfield=sapply(1:nrow(d),function(x){
  return(sqrt(sum((c(d$CN_seston[x],d$CP_seston[x],d$NP_seston[x])-
                     c(106/16,106,16))**2)))
})

write.table(d,"../Table/Simulation_random_parameters.csv",sep=";")
