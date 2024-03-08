packages = c(
  "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "igraph", "ggforce",
  "JuliaCall", "diffeqr", "phaseR", "ggtext", "viridis", "rootSolve",
  "ggquiver", "scales", "boot", "RColorBrewer", "ggnewscale","ggpattern","FME",
  "FactoMineR","factoextra"
)

# install pacakges if not installed already
install.packages(setdiff(packages, rownames(installed.packages())))


x = c(
  "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "igraph", "ggforce",
  "JuliaCall", "diffeqr", "phaseR", "ggtext", "viridis", "rootSolve",
  "ggquiver", "scales", "boot", "RColorBrewer", "ggnewscale","ggpattern","FME",
  "FactoMineR","factoextra"
)

lapply(x, require, character.only = TRUE)

# julia_setup()
# de = diffeq_setup()
# julia_library("DifferentialEquations")

the_theme = theme_classic() + theme(
  legend.position = "bottom",
  strip.background = element_rect(fill = "#CCE8D8"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 8),
  legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)



the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))


pal=colorRampPalette(c("#5C7ECC","#32C4E2","#7DCE9A","#B4CE7D","#EFDD35","#FFB700"))
pal_aqua=colorRampPalette((c("white","#CFDCDE","#A4DEE6","#49CADC","#3B7FE0","#193C82","#060D61")))
pal_terr=colorRampPalette((c("white","#D8ECCD","#BBE0A7","#86D45C","#3BA23B","#066F06","#033E03")))


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))



## Creating folders

dir.create("../Figures/", showWarnings = FALSE)
dir.create("../Figures/SI", showWarnings = FALSE)
dir.create("../Table/", showWarnings = FALSE)

`%!in%` = Negate(`%in%`)


# ------------------------------ Parameters and dynamics related function ---- 
# >> Parameters ----

eval_string = function(string, ...) eval(parse(text = string), ...)

Get_parameter_lake = function(kP, kN, muF, muO, dO, dF, beta_F, alpha_F, beta_O, alpha_O,        # plankton
                              eB, aN, aP, aD, dB, m, alpha_B, beta_B,                    # decomposers
                              IN, IP, ID, lN, lP, lD, beta_allo, alpha_allo, same_stoichio,    # allochtonous flows
                              functional_response_phyto, DC_,sO,sB,sF
) { 
  return(
    list(
      
      #Planktonic species
      
      kP      = kP,      # half-saturating constant P intake
      kN      = kN,      # half-saturating constant N intake
      muF      = muF,      # maximal growth rate fixers
      muO      = muO,      # maximal growth rate non-fixers
      dF      = dF,      # death rate fixers
      dO      = dO,      # death rate non-fixers
      beta_F  = beta_F,  # P:C fixers
      beta_O  = beta_O,  # P:C non-fixers
      alpha_F = alpha_F, # N:C fixers
      alpha_O = alpha_O, # N:C non-fixers
      
      #Decomposers and decomposition
      
      eB  = eB, # Decomposers efficiency
      aN  = aN, # Decomposers attack rate on nitrogen
      aP  = aP, # Decomposers attack rate on phosphorous
      aD  = aD, # Decomposers attack rate on detritus
      m  = m,   # Mineralization rate 
      dB  = dB, # Decomposers death rate
      alpha_B  = alpha_B, # Decomposers N:C ratio
      beta_B   = beta_B, # Decomposers N:C ratio
      
      # Allochtonous flows
      
      IN = IN, # Allochtonous inputs of nitrogen
      IP = IP, # Allochtonous inputs of phosphorous
      ID = ID, # Allochtonous inputs of detritus
      lN = lN, # Leaching rate of nitrogen
      lP = lP, # Leaching rate of phosphorous
      lD = lD, # Leaching rate of detritus
      beta_allo  = beta_allo,  # P:C ratio of allochtonous flows
      alpha_allo = alpha_allo,  # N:C ratio of allochtonous flows
      same_stoichio = same_stoichio, #Does allochtonous flows have the same stoichio as in the ecosystem
      functional_response_phyto=functional_response_phyto, # Type I or type II or DC
      DC_=DC_,
      sB = sB, # self-regulation of decomposers 
      sF = sF, # self-regulation of fixers
      sO = sO  # self-regulation of non-fixers
    )
  )
}

Get_classical_param_lake = function() {
  # We take 3 different scenarios
  `%!in%` = Negate(`%in%`)
  param = Get_parameter_lake(
    kP      = .1,
    kN      = .1,
    muF     = .24,
    muO     = .25,
    dF      = .2,
    dO      = .2,
    beta_F  = 1/(16*8),
    beta_O  = 1/(16*8),
    beta_B  = .15/10,
    alpha_F = 1/8,
    alpha_O = 1/8,
    alpha_B = .15,
    eB = .5,
    aN = .4,
    aP = .4,
    aD = .83,
    m  = .5,
    dB = .1,
    IN = 5,
    IP = 1,
    ID = 1,
    lN = 1,
    lP = 1,
    lD = 2,
    beta_allo  = .05/10,
    alpha_allo = .05,
    same_stoichio=F,
    functional_response_phyto=1,
    DC_=F,
    sF=1,
    sO=1,
    sB=1
  )
  
  return(param)
}

Get_initial_values = function(param,random_=F) {
  
  if (random_){
    state = c("BC" = runif(1,0,5), "BN" = 0, "BP" = 0, "FC" = runif(1,0,5),"FN" = 0,"FP" = 0, 
              "OC" = runif(1,0,5),"ON" = 0,"OP" = 0, "DC" = 2, "DN" = .1, "DP" = .01,
              "N" = 3, "P" = 2)
  }else{
    state = c("BC" = 4, "BN" = 0, "BP" = 0, "FC" = 2,"FN" = 0,"FP" = 0, 
              "OC" = 2,"ON" = 2,"OP" = 2, "DC" = 2, "DN" = .1, "DP" = .01,
              "N" = 3, "P" = 2)
  }
 
  state["BN"] = param$alpha_B * state["BC"]
  state["BP"] = param$beta_B  * state["BC"]
  state["FN"] = param$alpha_F * state["FC"]
  state["FP"] = param$beta_F  * state["FC"]
  state["ON"] = param$alpha_O * state["OC"]
  state["OP"] = param$beta_O  * state["OC"]
  return(state)
}


# >> Dynamics ----


# ode_lake_CNP_regulation = julia_eval("
# 
# function ode_lake_CNP_regulation(du, u, p, t)
# 
#     u[u .< 10^-8] .= 0
#     kP, kN, muF, muO, dF,dO, beta_F, beta_O,alpha_F, alpha_O,eB, aN, aP, aD, m,dB, alpha_B, beta_B, IN, IP, ID, lN, lP, lD, beta_allo, alpha_allo,same_stoichio,functional_response_phyto,DC_, sB, sF, sO = p
#     BC, BN, BP, FC,FN,FP, OC,ON,OP, DC, DN, DP, N, P  = u
# 
# 
#     if (DC_==1) #donnor controlled
#       gO_NP = copy((muO * min(P, N))) # growth rate non-fixer
#       gF_P  = copy((muF * P))             # growth rate fixer
#       uptake_D = copy(eB * aD * DC)             # uptake detritus
#       uptake_N = copy(aN * N)                   # uptake nitrogen
#       uptake_P = copy(aP * P)                   # uptake phosphorous
# 
#     else
#       if functional_response_phyto==1
#          gO_NP = copy(muO * min(P, N) * OC)  # growth rate non-fixer
#          gF_P  = copy(muF * P * FC)             # growth rate fixer
#       else #type 2
#          gO_NP = copy(muO * min(P/(kP+P),N/(kN+N)) * OC) # growth rate non-fixer
#          gF_P  = copy((muF * P)/(kP+P) * FC)             # growth rate fixer
#       end
#       uptake_D = copy(eB * aD * DC * BC)             # uptake detritus
#       uptake_N = copy(aN * N * BC)                   # uptake nitrogen
#       uptake_P = copy(aP * P * BC)                   # uptake phosphorous
# 
#     end
# 
#     beta_D  = copy(DP/DC)                     # P:C ratio detritus in the lake
#     alpha_D = copy(DN/DC)                     # N:C ratio detritus in the lake
# 
# 
#     # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
#     # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
# 
#     if same_stoichio == 1
#        alpha_allo = copy(alpha_D)
#        beta_allo  = copy(beta_D)
#     end
# 
#     phi_P_C_lim = copy((uptake_D * (beta_B - beta_D)) / beta_B)
#     phi_P_N_lim = copy(((beta_B - beta_D) * alpha_B / ( (alpha_B - alpha_D) * beta_B)) * uptake_N)
#     phi_P_P_lim = copy(uptake_P)
# 
#     phi_P =  copy(min.(phi_P_C_lim,phi_P_N_lim,phi_P_P_lim))      # immobilization/mineralization of P
#     phi_N =  copy((( (alpha_B - alpha_D) * beta_B) / ((beta_B - beta_D) * alpha_B)) * phi_P)      # immobilization/mineralization of N
#     phi_D =  copy(beta_B  * phi_P / (beta_B  - beta_D))      # decomposition of detritus
# 
#     #Now the dynamics :
# 
#     #DECOMPOSERS
#     du[1] = phi_D                                       - dB * BC           - m * BC  -           sB * BC * BC          #BC
#     du[2] = phi_D * alpha_D + phi_N * alpha_B - alpha_B * dB * BC - alpha_B * m * BC  - alpha_B * sB * BC * BC          #BN
#     du[3] = phi_D * beta_D  + phi_P * beta_B  - beta_B  * dB * BC - beta_B  * m * BC  - beta_B  * sB * BC * BC          #BP
# 
#     #PLANKTON
#     du[4] =            gF_P  - dF * FC    -           sF * FC * FC                     #FC
#     du[5] = alpha_F * (gF_P  - dF * FC)   - alpha_F * sF * FC * FC                     #FN
#     du[6] = beta_F  * (gF_P  - dF * FC)   - beta_F  * sF * FC * FC                     #FP
#     du[7] =            gO_NP - dO * OC    -           sO * OC * OC                     #OC
#     du[8] = alpha_O * (gO_NP - dO * OC)   - alpha_O * sO * OC * OC                     #ON
#     du[9] = beta_O  * (gO_NP - dO * OC)   - beta_O  * sO * OC * OC                     #OP
# 
#     #DETRITUS
#     du[10] = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D                #DC
# 
#     du[11] = ID * alpha_allo - lD * DN + alpha_B * dB * BC +                   #DN
#       alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D
# 
#     du[12] = ID * beta_allo  - lD * DP + beta_B  * dB * BC +                   #DP
#       beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D
# 
#     #NITROGEN
#     du[13] = IN - lN * N - alpha_O * gO_NP  -                              #N
#       alpha_B * phi_N + alpha_B * m * BC
# 
#     #PHOSPHOROUS
#     du[14] = IP - lP * P - beta_O * gO_NP   - beta_F  * gF_P  -        #P
#       beta_B * phi_P  + beta_B  * m * BC
# 
# end
# 
# ")

# for (k in 1:length(param)){assign(names(param)[k],param[[k]])}



ode_lake_CNP_R = function(t,y,param){ 
  
  y[y < 10^-8] = 0 # prevent numerical problems
  
  with(as.list(c(y, param)),{
    
    if (DC_){ #donnor controlled
      
      gO_NP = (muO * min(P, N))             # growth rate non-fixer
      gF_P  = (muF * P)                     # growth rate fixer
      uptake_D = (eB * aD * DC)             # uptake detritus
      uptake_N = (aN * N)                   # uptake nitrogen
      uptake_P = (aP * P)                   # uptake phosphorous
      
    } else{
      if (functional_response_phyto==1){
        
        gO_NP = muO * min(P, N) * OC     # growth rate non-fixer
        gF_P  = muF * P * FC             # growth rate fixer
        
      }else{ #type 2
        
        gO_NP = muO * min(P/(kP+P),N/(kN+N)) * OC # growth rate non-fixer
        gF_P  = (muF * P)/(kP+P) * FC             # growth rate fixer
      }
      uptake_D = eB * aD * DC * BC             # uptake detritus
      uptake_N = aN * N * BC                   # uptake nitrogen
      uptake_P = aP * P * BC                   # uptake phosphorous
    }
    
    beta_D  = DP/DC                     # P:C ratio detritus allocht
    alpha_D = DN/DC                     # N:C ratio detritus
    
    if (same_stoichio){
      alpha_allo = alpha_D
      beta_allo  = beta_D
    }
    
    # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
    # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
    
    #immobilization/mineralization P
    phi_P_C_lim = (uptake_D * (beta_B - beta_D)) / beta_B
    phi_P_N_lim = ((beta_B - beta_D) * alpha_B / ( (alpha_B - alpha_D) * beta_B)) * uptake_N
    phi_P_P_lim = uptake_P
    
    phi_P =  min(c(phi_P_C_lim,phi_P_N_lim,phi_P_P_lim))      # immobilization/mineralization of P
    phi_N =  (( (alpha_B - alpha_D) * beta_B) / ((beta_B - beta_D) * alpha_B)) * phi_P      # immobilization/mineralization of N
    phi_D =  beta_B  * phi_P / (beta_B  - beta_D)      # decomposition of detritus
    
    
    
    #Now the dynamics :
    
    du=rep(NA,14)
    
    #DECOMPOSERS     
    du[1] = phi_D                             -           dB * BC           - m * BC  -           sB * BC * BC  #BC
    du[2] = phi_D * alpha_D + phi_N * alpha_B - alpha_B * dB * BC - alpha_B * m * BC  - alpha_B * sB * BC * BC  #BN
    du[3] = phi_D * beta_D  + phi_P * beta_B  - beta_B  * dB * BC - beta_B  * m * BC  - beta_B  * sB * BC * BC  #BP
    
    #PLANKTON
    du[4] =            gF_P  - dF * FC    -            sF * FC * FC     #FC
    du[5] = alpha_F * (gF_P  - dF * FC)   - alpha_F  * sF * FC * FC     #FN
    du[6] = beta_F  * (gF_P  - dF * FC)   - beta_F   * sF * FC * FC     #FP
    du[7] =            gO_NP - dO * OC    -            sO * OC * OC     #OC
    du[8] = alpha_O * (gO_NP - dO * OC)   - alpha_O  * sO * OC * OC     #ON
    du[9] = beta_O  * (gO_NP - dO * OC)   - beta_O   * sO * OC * OC     #OP
    
    #DETRITUS 
    du[10] = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D                #DC  
    
    du[11] = ID * alpha_allo - lD * DN + alpha_B * dB * BC +                   #DN
      alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D                                    
    
    du[12] = ID * beta_allo  - lD * DP + beta_B  * dB * BC +                   #DP
      beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D                                    
    
    #NITROGEN 
    du[13] = IN - lN * N - alpha_O * gO_NP  - #- alpha_F * gF_P * FC       #N                           
      alpha_B * phi_N + alpha_B * m * BC
    
    #PHOSPHOROUS 
    du[14] = IP - lP * P - beta_O  * gO_NP   - beta_F  * gF_P  -        #P                      
      beta_B * phi_P  + beta_B  * m * BC
    
    #sum(du)
    #IP+ID+IN+ID * alpha_allo+ID * beta_allo-m * BC-lN * N-lP * P-lD * DC-lD * DN-lD * DP
    
    list(du)
    
  })
}

ode_lake_CNP_R_flexible = function(t,y,param){ 
  
  y[y < 10^-8] = 0 # prevent numerical problems
  
  with(as.list(c(y, param)),{
    
    if (DC_){ #donnor controlled
      
      gO_NP = (muO * min(P, N))             # growth rate non-fixer
      gF_P  = (muF * P)                     # growth rate fixer
      uptake_D = (eB * aD * DC)             # uptake detritus
      uptake_N = (aN * N)                   # uptake nitrogen
      uptake_P = (aP * P)                   # uptake phosphorous
      
    } else{
      if (functional_response_phyto==1){
        
        gO_NP = muO * min(P, N) * OC     # growth rate non-fixer
        gF_P  = muF * P * FC             # growth rate fixer
        
      }else{ #type 2
        
        gO_NP = muO * min(P/(kP+P),N/(kN+N)) * OC # growth rate non-fixer
        gF_P  = (muF * P)/(kP+P) * FC             # growth rate fixer
      }
      uptake_D = eB * aD * DC * BC             # uptake detritus
      uptake_N = aN * N * BC                   # uptake nitrogen
      uptake_P = aP * P * BC                   # uptake phosphorous
    }
    
    beta_D  = DP/DC                     # P:C ratio detritus allocht
    alpha_D = DN/DC                     # N:C ratio detritus
    
    if (same_stoichio){
      alpha_allo = alpha_D
      beta_allo  = beta_D
    }
    
    # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
    # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
    
    #immobilization/mineralization P
    phi_P_C_lim = (uptake_D * (beta_B - beta_D)) / beta_B
    phi_P_N_lim = ((beta_B - beta_D) * alpha_B / ( (alpha_B - alpha_D) * beta_B)) * uptake_N
    phi_P_P_lim = uptake_P
    
    phi_P =  min(c(phi_P_C_lim,phi_P_N_lim,phi_P_P_lim))      # immobilization/mineralization of P
    phi_N =  (( (alpha_B - alpha_D) * beta_B) / ((beta_B - beta_D) * alpha_B)) * phi_P      # immobilization/mineralization of N
    phi_D =  beta_B  * phi_P / (beta_B  - beta_D)      # decomposition of detritus
    
    
    
    #Now the dynamics :
    
    du=rep(NA,14)
    
    #DECOMPOSERS     
    du[1] = phi_D                             -           dB * BC           - m * BC  -           sB * BC * BC  #BC
    du[2] = phi_D * alpha_D + phi_N * alpha_B - alpha_B * dB * BC - alpha_B * m * BC  - alpha_B * sB * BC * BC  #BN
    du[3] = phi_D * beta_D  + phi_P * beta_B  - beta_B  * dB * BC - beta_B  * m * BC  - beta_B  * sB * BC * BC  #BP
    
    #PLANKTON
    du[4] =            gF_P  - dF * FC    -            sF * FC * FC     #FC
    du[5] = alpha_F * (gF_P  - dF * FC)   - alpha_F  * sF * FC * FC     #FN
    du[6] = beta_F  * (gF_P  - dF * FC)   - beta_F   * sF * FC * FC     #FP
    du[7] =            gO_NP - dO * OC    -            sO * OC * OC     #OC
    du[8] = alpha_O * (gO_NP - dO * OC)   - alpha_O  * sO * OC * OC     #ON
    du[9] = beta_O  * (gO_NP - dO * OC)   - beta_O   * sO * OC * OC     #OP
    
    #DETRITUS 
    du[10] = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D                #DC  
    
    du[11] = ID * alpha_allo - lD * DN + alpha_B * dB * BC +                   #DN
      alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D                                    
    
    du[12] = ID * beta_allo  - lD * DP + beta_B  * dB * BC +                   #DP
      beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D                                    
    
    #NITROGEN 
    du[13] = IN - lN * N - alpha_O * gO_NP  - #- alpha_F * gF_P * FC       #N                           
      alpha_B * phi_N + alpha_B * m * BC
    
    #PHOSPHOROUS 
    du[14] = IP - lP * P - beta_O  * gO_NP   - beta_F  * gF_P  -        #P                      
      beta_B * phi_P  + beta_B  * m * BC
    
    #sum(du)
    #IP+ID+IN+ID * alpha_allo+ID * beta_allo-m * BC-lN * N-lP * P-lD * DC-lD * DN-lD * DP
    
    list(du)
    
  })
}

Compute_ode = function(state,             # initial biomass and densities
                       param,             # parameter list
                       TRESH = 1e-5,      # threshold for extinction
                       n_time = 10000,    # maximal time 
                       julia=T,           # if F, the dynamics are computed on R
                       solver="lsoda",
                       model= "Tyrrell"
) {
  
  if (julia){
    
    #if julia, then all parameters and initial states are pushed in julia and the temporal dynamics are computed
    
    julia_assign("state", state)
    julia_assign("p", unlist(param))
    
    tspan = c(0, n_time) # to avoid long transient
    julia_assign("tspan", tspan)
    
    prob = julia_eval("ODEProblem(ode_lake_CNP_regulation, state, tspan, p)")
    
    sol = de$solve(prob, de$Tsit5())
    d = as.data.frame(t(sapply(sol$u, identity)))
    d$time = sol$t
    d=d[,c(ncol(d),1:(ncol(d)-1))]
  } else {
    
    times=seq(0,n_time,1)
    if (model=="Tyrrell"){
      d=as.data.frame(ode(y=state,times=times,func=ode_lake_CNP_R,
                          parms = param,method = solver))
    }else{
      param$Fix=.5
      d=as.data.frame(ode(y=state,times=times,func=ode_lake_CNP_R_Koffel,
                          parms = param,method = solver))
    }
    
  } 
  
  final_point = as.numeric(d[nrow(d), -1])
  final_point[which(final_point < TRESH)] = 0
  d[nrow(d), -1] = final_point
  colnames(d) = c("Time",
                  "Decomposers_C","Decomposers_N","Decomposers_P",
                  "Fixers_C","Fixers_N","Fixers_P","Non_fixers_C","Non_fixers_N","Non_fixers_P",
                  "Detritus_C","Detritus_N","Detritus_P",
                  "Nitrogen","Phosphorous")
  
  return(as_tibble(d))
}

plot_dynamics = function(data, log_ = T) {
  
  colors=c("Decomposers_C"="#FDC070",
           "Fixers_C"="#76BB70",
           "Non_fixers_C"="#5D9BF3",
           "Detritus_C"="brown",
           "NC_Detritus"="black",
           "PC_Detritus"="grey",
           "Nitrogen"="#4A50FF",
           "Phosphorous"="#926CAB")
  p = ggplot(data%>%
               add_column(., 
                          NC_Detritus=.$Detritus_N/.$Detritus_C,
                          PC_Detritus=.$Detritus_P/.$Detritus_C
               )%>%
               melt(., id.vars="Time")%>%
               filter(., variable %in% c("Decomposers_C","Fixers_C",
                                         "Non_fixers_C","Detritus_C",
                                         "NC_Detritus","PC_Detritus",
                                         "Nitrogen","Phosphorous"))) +
    geom_line(aes(x = Time, y = value, color = variable), lwd = 1) +
    labs(x = "Time", y = "Patch density", color = "") +
    scale_color_manual(values = colors) +
    the_theme +
    theme(legend.box = "vertical")
  
  if (log_) {
    p = p + scale_x_log10()+scale_y_log10()
  }
  
  return(p)
}

Test_convergence=function(data){
  
  data=melt(data[-which(data$Time<max(data$Time)-500)],id.vars=c("Time"))
  d_slope=tibble()
  for (name_var in unique(data$variable)){
    
    data_filtered=filter(data,variable==name_var)
    
    summary_fit=summary(lm(value~Time,data=data_filtered))
    
    d_slope=rbind(d_slope,tibble(Slope=summary_fit$coefficients[2,1],
                                 Variable=name_var,
                                 pvalue=summary_fit$coefficients[2,4]))
  }
  
  test_convergence=ifelse(all(abs(d_slope$Slope)<.05) & all(d_slope$pvalue>.05),"Converged","Not converged")
  
  return(list(Convergence_test=test_convergence,Slopes=d_slope))
}

# ------------------------------ Equilibrium related function ---- 

Get_equilibrium = function(data, param) {
  data=as.data.frame(data,as_tibble(t(unlist(param))))
  return(data[nrow(data),-1])
}

Get_limitation = function(data, param) {
  
  if (nrow(data)>1){data=Get_equilibrium(data,param)}
  
  alpha_D=data$Detritus_N / data$Detritus_C
  beta_D=data$Detritus_P / data$Detritus_C
  
  CN_threshold = (((param$alpha_B - alpha_D) / param$alpha_B) * param$eB * param$aD * data$Detritus_C) / (param$aN * data$Nitrogen)
  CP_threshold = (((param$beta_B - beta_D) / param$beta_B) * param$eB * param$aD * data$Detritus_C) / (param$aP * data$Phosphorous)
  NP_threshold = (param$aN * data$Nitrogen) / (param$aP * data$Phosphorous)
  
  if (CN_threshold<1 & CP_threshold<1){
    limit="C"
  }else if (CN_threshold<1 & CP_threshold>1){
    limit="P"
  }else if (CN_threshold>1 & CP_threshold<1){
    limit="N"
  }else if (CN_threshold>1 & CP_threshold>1 & NP_threshold<1){
    limit="N"
  }else{
    limit="P"
  }
  
  if (param$functional_response_phyto==2){
    NP_threshold_NF = (data$Nitrogen/ (data$Nitrogen+param$kN)) / (data$Phosphorous/ (data$Phosphorous+param$kP)) # *(param$beta_F/param$alpha_F) 
  }else{
    NP_threshold_NF = (data$Nitrogen) / (data$Phosphorous) # *(param$beta_F/param$alpha_F) 
  }
  
  limit_NF=ifelse(NP_threshold_NF<1,"N","P")
  
  return(list(NP_threshold = NP_threshold,
              CP_threshold = CP_threshold,
              CN_threshold = CN_threshold,
              Limitation_Decompo = limit,
              NP_threshold_NF=NP_threshold_NF,
              Limitation_NF=limit_NF))
}

Get_CNP = function(data, param) {
  
  if (nrow(data)>1){data=Get_equilibrium(data,param)}
  
  alpha_D=data$Detritus_N / data$Detritus_C
  beta_D=data$Detritus_P / data$Detritus_C
  
  NP_dissoluted = data$Nitrogen / data$Phosphorous
  
  #mean over organisms and detritus
  
  C_seston = data$Decomposers_C + data$Fixers_C + data$Non_fixers_C + data$Detritus_C
  N_seston = data$Decomposers_N + data$Fixers_N + data$Non_fixers_N + data$Detritus_N
  P_seston = data$Decomposers_P + data$Fixers_P + data$Non_fixers_P + data$Detritus_P
  
  CN_seston = C_seston/N_seston
  CP_seston = C_seston/P_seston
  NP_seston = N_seston/P_seston
  
  #fraction of decomposers over plankton species
  
  frac_decompo = data$Decomposers_C / (data$Decomposers_C + data$Fixers_C + data$Non_fixers_C)
  
  return(list(CN_seston = CN_seston,
              CP_seston = CP_seston,
              NP_seston = NP_seston,
              C_seston  = C_seston, 
              N_seston  = N_seston, 
              P_seston  = P_seston, 
              NP_dissoluted = NP_dissoluted,
              NC_detritus = alpha_D,
              PC_detritus = beta_D,
              Frac_decomp = frac_decompo))
}

Plot_with_limitation=function(data,driver="ID",y_variable){
  
  data=melt(data,measure.vars = driver,variable.name = "variable_driver",
            value.name = "value_driver")%>%
    melt(., measure.vars = y_variable,variable.name = "variable_y",
         value.name = "value_y")
  
  for (k in unique(data$Limitation_Decompo)){
    
    assign(paste0("min_",k),min(data$value_driver[which(data$Limitation_Decompo==k)]))
    assign(paste0("max_",k),max(data$value_driver[which(data$Limitation_Decompo==k)]))
    
  } 
  
  if (length(unique(data$Limitation_NF))==1){
    min_NF=min(data$value_driver)
    max_NF=max(data$value_driver)
  }else{
    for (k in unique(data$Limitation_NF)){
      
      assign(paste0("min_NF_",k),min(data$value_driver[which(data$Limitation_NF==k)]))
      assign(paste0("max_NF_",k),max(data$value_driver[which(data$Limitation_NF==k)]))
      
    } 
  }
  
  min_y=min(data$value_y)
  max_y=max(data$value_y)

  diff_y=max_y-min_y
  min_y_pattern=max_y+.025*diff_y
  max_y_pattern=max_y+.125*diff_y
  
  if ("P" %!in% unique(data$Limitation_Decompo) &
      "N" %!in% unique(data$Limitation_Decompo)){
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_C,xmax=max_C,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#D2B96F")+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }else{
      
      d_rect=tibble(min_x=c(min_NF_N,min_NF_P),
                    max_x=c(max_NF_N,max_NF_P),
                    min_y=min_y,
                    max_y=max_y,
                    pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_C,xmax=max_C,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#D2B96F")+
        geom_rect_pattern(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
    
  }else if ("P" %!in% unique(data$Limitation_Decompo) &
            "C" %!in% unique(data$Limitation_Decompo)){
    
    min_y
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_N,xmax=max_N,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#8EBAEF")+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }else{
      
      d_rect=tibble(min_x=c(min_NF_N,min_NF_P),
                    max_x=c(max_NF_N,max_NF_P),
                    min_y=min_y,
                    max_y=max_y,
                    pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_N,xmax=max_N,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#8EBAEF")+
        geom_rect_pattern(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
    
  }else if ("N" %!in% unique(data$Limitation_Decompo) &
            "C" %!in% unique(data$Limitation_Decompo)){
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_P,xmax=max_P,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#D7B2F9")+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }else{
      
      d_rect=tibble(min_x=c(min_NF_N,min_NF_P),
                    max_x=c(max_NF_N,max_NF_P),
                    min_y=min_y,
                    max_y=max_y,
                    pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=NULL,aes(xmin=min_P,xmax=max_P,ymin=min_y,ymax=max_y),
                  alpha=.3,fill="#D7B2F9")+
        geom_rect_pattern(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                          alpha=.3,fill="grey")+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    

  }else if ("P" %in% unique(data$Limitation_Decompo) & 
            "C" %in% unique(data$Limitation_Decompo) & 
            "N" %in% unique(data$Limitation_Decompo) ){
    
    
    color_graph=c("C"="#D2B96F","N"="#8EBAEF","P"="#D7B2F9")
    
    d_rect=tibble(min_x=c(min_C,min_N,min_P),
                  max_x=c(max_C,max_N,max_P),
                  min_y=min_y,
                  max_y=max_y,
                  color=c("C","N","P"))%>%
      dplyr::arrange(., min_x)
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                          alpha=.3)+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                          alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe"))

    }else{
      
      d_rect2=tibble(min_x=c(min_NF_N,min_NF_P),
                    max_x=c(max_NF_N,max_NF_P),
                    min_y=min_y,
                    max_y=max_y,
                    pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=d_rect2,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                                           alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
    
    
  }else if ("P" %in% unique(data$Limitation_Decompo) & 
            "C" %!in% unique(data$Limitation_Decompo) & 
            "N" %in% unique(data$Limitation_Decompo) ){
    
    
    color_graph=c("N"="#8EBAEF","P"="#D7B2F9")
    
    d_rect=tibble(min_x=c(min_N,min_P),
                  max_x=c(max_N,max_P),
                  min_y=min_y,
                  max_y=max_y,
                  color=c("N","P"))%>%
      dplyr::arrange(., min_x)
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                                        alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe"))
      
    }else{
      
      d_rect2=tibble(min_x=c(min_NF_N,min_NF_P),
                     max_x=c(max_NF_N,max_NF_P),
                     min_y=min_y,
                     max_y=max_y,
                     pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=d_rect2,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                                           alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
  }else if ("P" %in% unique(data$Limitation_Decompo) & 
            "C" %in% unique(data$Limitation_Decompo) & 
            "N" %!in% unique(data$Limitation_Decompo) ){
    
    
    color_graph=c("C"="#D2B96F","P"="#D7B2F9")
    
    d_rect=tibble(min_x=c(min_C,min_P),
                  max_x=c(max_C,max_P),
                  min_y=min_y,
                  max_y=max_y,
                  color=c("C","P"))%>%
      dplyr::arrange(., min_x)
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                                        alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe"))
      
    }else{
      
      d_rect2=tibble(min_x=c(min_NF_N,min_NF_P),
                     max_x=c(max_NF_N,max_NF_P),
                     min_y=min_y,
                     max_y=max_y,
                     pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=d_rect2,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                                           alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
  }else if ("P" %!in% unique(data$Limitation_Decompo) & 
            "C" %in% unique(data$Limitation_Decompo) & 
            "N" %in% unique(data$Limitation_Decompo) ){
    
    
    color_graph=c("C"="#D2B96F","N"="#8EBAEF")
    
    d_rect=tibble(min_x=c(min_C,min_N),
                  max_x=c(max_C,max_N),
                  min_y=min_y,
                  max_y=max_y,
                  color=c("C","N"))%>%
      dplyr::arrange(., min_x)
    
    if (length(unique(data$Limitation_NF))==1){
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=NULL,aes(xmin=min_NF,xmax=max_NF,ymin=min_y_pattern,ymax=max_y_pattern,pattern=data$Limitation_NF[1]),
                                        alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe"))
      
    }else{
      
      d_rect2=tibble(min_x=c(min_NF_N,min_NF_P),
                     max_x=c(max_NF_N,max_NF_P),
                     min_y=min_y,
                     max_y=max_y,
                     pattern=c("N","P"))%>%
        dplyr::arrange(., min_x)
      
      
      p=ggplot(NULL)+
        geom_rect(data=d_rect,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,fill=color),
                  alpha=.3)+
        geom_rect_pattern(data=d_rect2,aes(xmin=min_x,xmax=max_x,ymin=min_y_pattern,ymax=max_y_pattern,pattern=pattern),
                                           alpha=.3,fill="grey")+
        scale_fill_manual(values=color_graph)+
        scale_pattern_manual(values = c("N" = "none", "P" = "stripe")) 
      
    }
    
  }
  
  return(list(p=p+guides(fill="none"),data=data))
}

Get_phi_functions=function(data,param){
  
  if (nrow(data)>1){data=Get_equilibrium(data,param)}
  
  colnames(data)=c("BC","BN","BP","FC","FN","FP","OC","ON","OP",
                   "DC","DN","DP","N","P")
  
  for (k in 1:length(param)){assign(names(param)[k],param[[k]])}
  
  if (DC_){ #donnor controlled
    
    gO_NP = (muO * min(data$P, data$N))             # growth rate non-fixer
    gF_P  = (muF * data$P)                     # growth rate fixer
    uptake_D = (eB * aD * data$DC)             # uptake detritus
    uptake_N = (aN * data$N)                   # uptake nitrogen
    uptake_P = (aP * data$P)                   # uptake phosphorous
    
  } else{
    if (functional_response_phyto==1){
      
      gO_NP = muO * min(data$P, data$N) * data$OC     # growth rate non-fixer
      gF_P  = muF * data$P * data$FC             # growth rate fixer
      
    }else{ #type 2
      
      gO_NP = muO * min(data$P/(kP+data$P),data$N/(kN+data$N)) * data$OC # growth rate non-fixer
      gF_P  = (muF * data$P)/(kP+data$P) * data$FC             # growth rate fixer
    }
    uptake_D = eB * aD * data$DC * data$BC             # uptake detritus
    uptake_N = aN * data$N * data$BC                   # uptake nitrogen
    uptake_P = aP * data$P * data$BC                   # uptake phosphorous
  }
  
  beta_D  = data$DP/data$DC                     # P:C ratio detritus allocht
  alpha_D = data$DN/data$DC                     # N:C ratio detritus
  
  if (same_stoichio){
    alpha_allo = alpha_D
    beta_allo  = beta_D
  }
  
  # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
  # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
  
  #immobilization/mineralization P
  phi_P_C_lim = (uptake_D * (beta_B - beta_D)) / beta_B
  phi_P_N_lim = ((beta_B - beta_D) * alpha_B / ( (alpha_B - alpha_D) * beta_B)) * uptake_N
  phi_P_P_lim = uptake_P
  
  phi_P =  min(c(phi_P_C_lim,phi_P_N_lim,phi_P_P_lim))      # immobilization/mineralization of P
  phi_N =  (( (alpha_B - alpha_D) * beta_B) / ((beta_B - beta_D) * alpha_B)) * phi_P      # immobilization/mineralization of N
  phi_D =  beta_B  * phi_P / (beta_B  - beta_D)      # decomposition of detritus
  
  return(tibble(phi_P=phi_P,
                phi_N=phi_N,
                phi_D=phi_D))
}

Get_states=function(data){
  
  states=sapply(1:nrow(data),function(x){
    
    if (data$Decomposers_C[x] !=0 & data$Fixers_C[x] !=0 & data$Non_fixers_C[x] !=0){
      return("Coexistence")
    } else if (data$Decomposers_C[x] ==0 & data$Fixers_C[x] !=0 & data$Non_fixers_C[x] !=0){
      return("Coexistence F-NF")
    } else if (data$Decomposers_C[x] !=0 & data$Fixers_C[x] ==0 & data$Non_fixers_C[x] !=0){
      return("Coexistence D-NF")
    } else if (data$Decomposers_C[x] !=0 & data$Fixers_C[x] !=0 & data$Non_fixers_C[x] ==0){
      return("Coexistence D-F")
    } else if (data$Decomposers_C[x] ==0 & data$Fixers_C[x] ==0 & data$Non_fixers_C[x] !=0){
      return("F only")
    } else if (data$Decomposers_C[x] ==0 & data$Fixers_C[x] ==0 & data$Non_fixers_C[x] !=0){
      return("NF only")
    } else if (data$Decomposers_C[x] !=0 & data$Fixers_C[x] ==0 & data$Non_fixers_C[x] ==0){
      return("B only")
    } else{
      return("Extinction")
    }
  })
  
  return(states)
}

Add_CNP=function(d,param){
  
  d[,c("CN_seston","CP_seston","NP_seston")]=NA
  
  for (x in 1:nrow(d)){
    CNP_x=Get_CNP(d[x,],param)
    d$CN_seston[x]=CNP_x$CN_seston
    d$CP_seston[x]=CNP_x$CP_seston
    d$NP_seston[x]=CNP_x$NP_seston
  }
  
  return(d)
}

Plot_D=function(d){
  
  if ("CP_seston" %!in% colnames(d)){
    d=Add_CNP(d,param)
  }
  
  par(mfrow=c(3,3),mar=rep(2,4))
  plot(d$Non_fixers_C)
  plot(d$Fixers_C)
  plot(d$Decomposers_C)
  plot(d$Nitrogen)
  plot(d$Phosphorous)
  plot(d$Detritus_C)
  plot(d$CN_seston)
  plot(d$CP_seston)
  plot(d$NP_seston)
  par(mfrow=c(1,1))
}

# ------------------------------ Indirect effects ----

Indirect_effects=function(data,param){
  
  d_effect=tibble()
  
  if (nrow(data)>1){
    data=Get_equilibrium(data,param)
  }
  data=as.numeric(data);names(data)=c("BC","BN","BP","FC","FN",
                                      "FP","OC","ON","OP","DC","DN","DP","N","P")
  
  jacobian=rootSolve::jacobian.full(data,ode_lake_CNP_R,parms = param,pert = 1e-9)
  jacobian=jacobian[-which(colnames(jacobian) %in% c("BN","BP","FN","FP","ON","OP")),
                    -which(colnames(jacobian) %in% c("BN","BP","FN","FP","ON","OP"))]

  colnames(jacobian)=rownames(jacobian)=c("B","F","O","DC","DN","DP","N","P")

  
  if (any(rowSums(jacobian)==0 | colSums(jacobian)==0)){
    to_remove_row=which(rowSums(jacobian)==0)
    to_remove_col=which(colSums(jacobian)==0)
    jacobian=jacobian[-to_remove_row,-to_remove_col]
    
  }
    
  sensitivity=-solve(jacobian)
  
  save_J=jacobian
  save_S=sensitivity
  
  #inflow to abundance
  diag(sensitivity)=NA
  
  for (name_organism1 in c("F","B","O")){
    for (name_organism2 in c("F","B","O")[-which(c("F","B","O")==name_organism1)]){
      if (name_organism1 %in% colnames(sensitivity) & 
          name_organism2 %in% colnames(sensitivity)){
        assign(paste0("indirect_",name_organism1,"_on_",name_organism2),
               sensitivity[name_organism2,name_organism1]
               )
      }else{
        assign(paste0("indirect_",name_organism1,"_on_",name_organism2),
               NA
        )
        
      }
    }
  }
  

  #Since there is no direct effects between species, net effects correspond to indirect effects
  
  d_effect=rbind(d_effect,tibble(indirect_B_on_F=indirect_B_on_F,
                                 indirect_F_on_B=indirect_F_on_B,
                                 indirect_O_on_F=indirect_O_on_F,
                                 indirect_F_on_O=indirect_F_on_O,
                                 indirect_B_on_O=indirect_B_on_O,
                                 indirect_O_on_B=indirect_O_on_B,
                                 Type="I_to_A"
  ))  
  
  
  #abundance to inflow
  sensitivity=save_S
  for(i in 1:nrow(jacobian)){
    for(j in 1:ncol(jacobian)){
      sensitivity[i,j]=save_S[i,j]/(save_S[i,i]*save_S[j,j]-save_S[i,j]*save_S[j,i])
    }}
  diag(sensitivity)=NA
  
  
  for (name_organism1 in c("F","B","O")){
    for (name_organism2 in c("F","B","O")[-which(c("F","B","O")==name_organism1)]){
      if (name_organism1 %in% colnames(sensitivity) & 
          name_organism2 %in% colnames(sensitivity)){
        assign(paste0("indirect_",name_organism1,"_on_",name_organism2),
               sensitivity[name_organism2,name_organism1]
        )
      }else{
        assign(paste0("indirect_",name_organism1,"_on_",name_organism2),
               NA
        )
        
      }
    }
  }
  
  d_effect=rbind(d_effect,tibble(indirect_B_on_F=indirect_B_on_F,
                                 indirect_F_on_B=indirect_F_on_B,
                                 indirect_O_on_F=indirect_O_on_F,
                                 indirect_F_on_O=indirect_F_on_O,
                                 indirect_B_on_O=indirect_B_on_O,
                                 indirect_O_on_B=indirect_O_on_B,
                                 Type="A_to_I"
  ))  
  
  return(d_effect)
}

Plot_net_effects=function(data,param,type="A_to_I",compute_effects=T){
  
  if (compute_effects){
    data=cbind(Indirect_effects(data,param)%>%
      filter(., Type==type),
      Get_CNP(data,param),
      Get_limitation(data,param)
      )
  }
  
  Matrix_indirect=matrix(c(0,data$indirect_B_on_O,data$indirect_B_on_F,
                           data$indirect_O_on_B,0,data$indirect_O_on_F,
                           data$indirect_F_on_B,data$indirect_F_on_O,0
                           ),
                         3,3,byrow = T)
  colnames(Matrix_indirect)=rownames(Matrix_indirect)=c("B","NF","F")
  
  graph_indirect=graph_from_adjacency_matrix(Matrix_indirect,weighted = T,mode = "directed",diag = F)  
  
  E(graph_indirect)$color=ifelse(E(graph_indirect)$weight<0,"#FF5B5B","#6880D0")
  E(graph_indirect)$size=E(graph_indirect)$weight
  V(graph_indirect)$color=c(ifelse(data$Limitation_Decompo=="N","#B9D6FF",
                                   ifelse(data$Limitation_Decompo=="P","#F0B8FF","#EACDA0")),
                            ifelse(data$Limitation_NF=="N","#B9D6FF","#F0B8FF"),
                            "#F0B8FF")
  
  plot(graph_indirect,layout=layout.circle,edge.curved=.3,
       vertex.size=55,edge.width=10*E(graph_indirect)$weight,
       label.cex=25)
}
