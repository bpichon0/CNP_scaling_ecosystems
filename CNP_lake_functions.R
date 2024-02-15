packages = c(
    "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "igraph", "ggforce",
    "JuliaCall", "diffeqr", "phaseR", "ggtext", "viridis", "rootSolve",
    "ggquiver", "scales", "boot", "RColorBrewer", "ggnewscale"
)

# install pacakges if not installed already
install.packages(setdiff(packages, rownames(installed.packages())))


x = c(
    "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", "igraph", "ggforce",
    "JuliaCall", "diffeqr", "phaseR", "ggtext", "viridis", "rootSolve",
    "ggquiver", "scales", "boot", "RColorBrewer", "ggnewscale"
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


## Creating folders

dir.create("../Figures/", showWarnings = FALSE)
dir.create("../Figures/SI", showWarnings = FALSE)
dir.create("../Table/", showWarnings = FALSE)




# ------------------------------ Parameters and dynamics related function ---- 
# >> Parameters ----

eval_string = function(string, ...) eval(parse(text = string), ...)


Get_parameter_lake = function(kP, kN, µF, µO, dO, dF, beta_F, alpha_F, beta_O, alpha_O,        # plankton
                               eB, aN, aP, aD, dB, m, alpha_B, beta_B,                    # decomposers
                               IN, IP, ID, lN, lP, lD, beta_allo, alpha_allo, same_stoichio    # allochtonous flows
) { 
  return(
    list(
      
      #Planktonic species
      
      kP      = kP,      # half-saturating constant P intake
      kN      = kN,      # half-saturating constant N intake
      µF      = µF,      # maximal growth rate fixers
      µO      = µO,      # maximal growth rate non-fixers
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
      same_stoichio = same_stoichio #Does allochtonous flows have the same stoichio as in the ecosystem
    )
  )
}


Get_classical_param_lake = function() {

    # We take 3 different scenarios
    `%!in%` = Negate(`%in%`)
    param = Get_parameter_lake(
      kP      = .1,
      kN      = .1,
      µF      = 1,
      µO      = 1,
      dF      = .1,
      dO      = .1,
      beta_F  = .05/16,
      beta_O  = .05/16,
      beta_B  = .15/10,
      alpha_F = .05,
      alpha_O = .05,
      alpha_B = .15,
      eB = .5,
      aN = .1,
      aP = .1,
      aD = .83,
      m  = .5,
      dB = .1,
      IN = 1,
      IP = 1,
      ID = 1,
      lN = .1,
      lP = .1,
      lD = .1,
      beta_allo  = .05/20,
      alpha_allo = .05,
      same_stoichio=F
    )
    
    return(param)
}



Get_initial_values = function(param) {
  state = c("BC" = 2, "BN" = 0, "BP" = 0, "FC" = 2, "OC" = 2, "DC" = 2, "DN" = .1, "DP" = .1, "N" = 2, "P" = 2)
  state["BN"] = param$alpha_B * state["BC"]
  state["BP"] = param$beta_B  * state["BC"]
  return(state)
}


# >> Dynamics ----
# 
# ode_lake_CNP = julia_eval("
# 
# function ode_lake_CNP(du, u, p, t)
# 
#     u[u .< 1e-5] .= 0
#     kP, kN, µF, µO, dO, dF, beta_F, alpha_F, beta_O, alpha_O,eB, aN, aP, aD, dB, m, alpha_B, beta_B, IN, IP, ID, lN, lP, lD, beta_allo, alpha_allo,same_stoichio = p
#     BC, BN, BP, FC, OC, DC, DN, DP, N, P  = u
#     
#     gO_NP = copy(µO * min(P/(kP+P),N/(kN+N))) # growth rate non-fixer
#     gF_P  = copy((µF * P)/(kP+P))             # growth rate fixer
#     
#     uptake_D = copy(eB * aD * DC)             # uptake detritus
#     uptake_N = copy(aN * N)                   # uptake nitrogen
#     uptake_P = copy(aP * P)                   # uptake phosphorous
#     
#     beta_D  = copy(DP/DC)                      # P:C ratio detritus in the lake
#     alpha_D = copy(DN/DC)                      # N:C ratio detritus in the lake
#     
#     
#     # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
#     # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
#     
#     if same_stoichio
#        alpha_allo = copy(alpha_D)
#        beta_allo  = copy(beta_D)
#     end
#     #immobilization/mineralization P under the different limitations
#     
#     phi_P_C = copy((uptake_D * (beta_B - beta_D)) / beta_B)
#     phi_P_N = copy(((beta_B - beta_D) * alpha_D / ( (alpha_B - alpha_D) * beta_B)) * uptake_N)
#     phi_P_P = copy(uptake_P)
#     
#     
#     
#     #immobilization/mineralization N under the different limitations
#     
#     phi_N_C = copy((uptake_D * (beta_B - alpha_D)) / alpha_B)
#     phi_N_N = copy(uptake_N)
#     phi_N_P = copy((( (alpha_B - alpha_D) * beta_B) / (beta_B - beta_D) * alpha_D) * uptake_P)
#     
#     
#     
#     #decomposition C under the different limitations
#     
#     phi_D_C = copy(uptake_D)
#     phi_D_N = copy(alpha_B * uptake_P / (alpha_B - alpha_D))
#     phi_D_P = copy(beta_B  * uptake_P / (beta_B  - beta_D))
#     
#     
#     phi_P = copy( min(phi_P_C,phi_P_N,phi_P_P) )     # immobilization/mineralization of P
#     phi_N = copy( min(phi_N_C,phi_N_N,phi_N_P) )     # immobilization/mineralization of N
#     phi_D = copy( min(phi_D_C,phi_D_N,phi_D_P) )     # decomposition of detritus
# 
# 
#     #Now the dynamics :
#     
#     #DECOMPOSERS     
#     du[1] = phi_D                              - dB * BC - m * BC             #BC
#     du[2] = phi_D * alpha_D + alpha_B * (phi_N - dB * BC - m * BC)            #BN
#     du[3] = phi_D * beta_D  + beta_B  * (phi_P - dB * BC - m * BC)            #BP
#     
#     #PLANKTON
#     du[4] = gF_P  * FC - dF * FC                                              #FC
#     du[5] = gO_NP * OC - dO * OC                                              #OC
# 
#     #DETRITUS 
#     du[6]  = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D               #DC  
#     
#     du[7] = ID * alpha_allo - lD * DN + alpha_B * dB * BC +                   #DN
#     alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D                                    
#     
#     du[8] = ID * beta_allo - lD * DP + beta_B * dB * BC +                     #DP
#     beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D                                    
# 
#     #NITROGEN 
#     du[9] = IN - lN * N + alpha_O * gO_NP * OC - alpha_F * gF_P * FC -        #N                           
#              alpha_B * phi_N + alpha_B * m * BC
# 
#     #PHOSPHOROUS 
#     du[10] = IP - lP * P + beta_O * gO_NP * OC - beta_F * gF_P * FC -         #P                      
#              beta_B * phi_P + beta_B * m * BC
# 
# end
# 
# ")

ode_lake_CNP_R = function(t,y,param){ 
  
  with(as.list(c(y, param)),{
    
    gO_NP = µO * min(P/(kP+P),N/(kN+N)) # growth rate non-fixer
    gF_P  = (µF * P)/(kP+P)             # growth rate fixer
    
    uptake_D = eB * aD * DC * BC        # uptake detritus
    uptake_N = aN * N * BC              # uptake nitrogen
    uptake_P = aP * P * BC              # uptake phosphorous
    
    beta_D  = DP/DC                     # P:C ratio detritus allocht
    alpha_D = DN/DC                     # N:C ratio detritus
    
    if (same_stoichio){
      alpha_allo = alpha_D
      beta_allo  = beta_D
    }
    
    # we below write phi_i_j the decomposer function (immobilization/mineralization or decomposition) under the limitation j
    # phi_P_C would for instance be the immobilization/mineralization of P under C-limitation of decomposers
    
    #immobilization/mineralization P
    phi_P_C = (uptake_D * (beta_B - beta_D)) / beta_B
    phi_P_N = ((beta_B - beta_D) * alpha_D / ( (alpha_B - alpha_D) * beta_B)) * uptake_N
    phi_P_P = uptake_P
    
    #immobilization/mineralization N
    phi_N_C = (uptake_D * (beta_B - alpha_D)) / alpha_B
    phi_N_N = uptake_N
    phi_N_P = (( (alpha_B - alpha_D) * beta_B) / (beta_B - beta_D) * alpha_D) * uptake_P
    
    #decomposition C
    phi_D_C = uptake_D
    phi_D_N = alpha_B * uptake_P / (alpha_B - alpha_D)
    phi_D_P = beta_B  * uptake_P / (beta_B  - beta_D)
    
    
    phi_P =  min(phi_P_C,phi_P_N,phi_P_P)      # immobilization/mineralization of P
    phi_N =  min(phi_N_C,phi_N_N,phi_N_P)      # immobilization/mineralization of N
    phi_D =  min(phi_D_C,phi_D_N,phi_D_P)      # decomposition of detritus
    
    #Now the dynamics :
    
    du=rep(NA,10)
    
    #DECOMPOSERS     
    du[1] = phi_D                              - dB * BC - m * BC             #BC
    du[2] = phi_D * alpha_D + alpha_B * (phi_N - dB * BC - m * BC)            #BN
    du[3] = phi_D * beta_D  + beta_B  * (phi_P - dB * BC - m * BC)            #BP
    
    #PLANKTON
    du[4] = gF_P  * FC - dF * FC                                              #FC
    du[5] = gO_NP * OC - dO * OC                                              #OC
    
    #DETRITUS 
    du[6] = ID - lD * DC + dB * BC + dF * FC + dO * OC - phi_D                #DC  
    
    du[7] = ID * alpha_allo - lD * DN + alpha_B * dB * BC +                   #DN
      alpha_F * dF * FC + alpha_O * dO * OC - phi_D * alpha_D                                    
    
    du[8] = ID * beta_allo  - lD * DP + beta_B  * dB * BC +                   #DP
      beta_F  * dF * FC + beta_O  * dO * OC - phi_D * beta_D                                    
    
    #NITROGEN 
    du[9] = IN - lN * N + alpha_O * gO_NP * OC - alpha_F * gF_P * FC -        #N                           
      alpha_B * phi_N + alpha_B * m * BC
    
    #PHOSPHOROUS 
    du[10] = IP - lP * P + beta_O * gO_NP * OC - beta_F  * gF_P * FC -        #P                      
      beta_B * phi_P  + beta_B  * m * BC
    
    list(du)
    
  })
}



Compute_ode = function(state,           # initial biomass and densities
                       param,           # parameter list
                       TRESH = 1e-5,    # threshold for extinction
                       n_time = 10000,  # maximal time 
                       julia=T,         # if F, the dynamics are computed on R
                       solver="lsoda"   # ode solver
) {
  
  if (julia){
    
    #if julia, then all parameters and initial states are pushed in julia and the temporal dynamics are computed
    
    julia_assign("state", state)
    julia_assign("p", unlist(param))
    
    tspan = c(0, n_time) # to avoid long transient
    julia_assign("tspan", tspan)
    
    prob = julia_eval("ODEProblem(ode_lake_CNP, state, tspan, p)")
    
    sol = de$solve(prob, de$Tsit5())
    d = as.data.frame(t(sapply(sol$u, identity)))
    d$time = sol$t
    
  } else {
    
    times=seq(0,n_time,1)
    d=as.data.frame(ode(y=state,times=times,func=ode_lake_CNP_R,
                               parms = param,method = solver))
  } 
  
  final_point = as.numeric(d[nrow(d), -1])
  final_point[which(final_point < TRESH)] = 0
  d[nrow(d), -1] = final_point
  colnames(d) = c("Time",
                  "Decomposers_C","Decomposers_N","Decomposers_P",
                  "Fixers_C","Non_fixers_C",
                  "Detritus_C","Detritus_N","Detritus_P",
                  "Nitrogen","Phosphorous")
  
  return(as_tibble(d))
}



plot_dynamics = function(data, log_ = T) {
    colors = c("Consumers" = "darkorange", "Producers" = "green3", "Nitrogen" = "darkorchid2", "Detritus" = "brown", "Top predators" = "blue")

    the_theme = theme_classic() + theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "#CCE8D8"),
        strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
        strip.text.x = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10)
    )

    data = gather(data, variable, value, -Time)
    data$foodweb = Get_foodweb(data$variable)
    data$trophic_level = Get_trophic_level(data$variable)
    data$Resources = Get_type_resources(data$variable)
    p = ggplot(data) +
        geom_line(aes(x = Time, y = value, color = trophic_level, linetype = Resources), lwd = 1) +
        ylim(0, max(data$value)) +
        labs(x = "Time", y = "Patch density", color = "Trophic level", linetype = "Resource") +
        scale_color_manual(values = colors) +
        facet_grid(. ~ foodweb) +
        the_theme +
        theme(legend.box = "vertical")

    if (log_) {
        p = p + scale_x_log10()
    }

    return(p)
}



# ------------------------------ Equilibrium related function ---- 

Get_equilibrium = function(data, param, consumers = F) {
    n_begin = ifelse(consumers, 19, 15)

    data_mean = as_tibble(t(colMeans(data[(nrow(data) - 1000):nrow(data), -1])))
    data_with_param = cbind(data_mean, matrix(unlist(param), ncol = length(param), nrow = 1))
    colnames(data_with_param)[n_begin:ncol(data_with_param)] = names(param)

    return(list(Eq = data_with_param))
}

Get_type_resources = function(vector) {
    return(unlist(lapply(vector, function(x) {
        if (x %in% c("Top_cons_H_N", "Top_cons_C_N", "Herbivores_N", "Consumers_N", "Decomposers_N", "Plants_N", "Nitrogen_T_N", "Nitrogen_A_N", "Detritus_T_N", "Detritus_A_N")) {
            return("N")
        } else {
            return("C")
        }
    })))
}

Get_limitation = function(data, param) {
    if (class(data) == "numeric") {
        data = as_tibble(t(data))
        colnames(data) = c("Herbivores_C", "Herbivores_N", "Consumers_C", "Consumers_N", "Plants_C", "Plants_N", "Decomposers_C", "Decomposers_N", "Detritus_T_C", "Detritus_T_N", "Detritus_A_C", "Detritus_A_N", "Nitrogen_T_N", "Nitrogen_A_N")
    }

    if (((((param$rB - data$Detritus_A_N / data$Detritus_A_C) / param$rB) * param$eB * param$aBD * data$Detritus_A_C) / (param$aBN * data$Nitrogen_A_N)) < 1) {
        limit = "C-limited"
    }
    if (((((param$rB - data$Detritus_A_N / data$Detritus_A_C) / param$rB) * param$eB * param$aBD * data$Detritus_A_C) / (param$aBN * data$Nitrogen_A_N)) > 1) {
        limit = "N-limited"
    }

    ratio_C_N = ((((param$rB - data$Detritus_A_N / data$Detritus_A_C) / param$rB) * param$eB * param$aBD * data$Detritus_A_C) / (param$aBN * data$Nitrogen_A_N))

    return(list(Ratio = ratio_C_N, Limitation = limit))
}

Limitation_data = function(data) {
    data$Limitation = sapply(1:nrow(data), function(x) {
        param = data[x, 15:ncol(data)]

        Get_limitation(data[x, -1], param)$Ratio
    })
    return(data)
}

State_at_equilibrium = function(data, consumers = F) {
    data = round(data, 5)
    if (consumers == F) {
        if (data$Plants_C == 0 & data$Decomposers_C != 0) {
            state = "Terrestrial_extinct"
        } else if (data$Plants_C != 0 & data$Decomposers_C == 0) {
            state = "Aquatic_extinct"
        } else if (data$Plants_C != 0 & data$Herbivores_C != 0 & data$Decomposers_C != 0 & data$Consumers_C != 0) {
            state = "Coexistence"
        } else if (data$Plants_C != 0 & data$Herbivores_C == 0 & data$Decomposers_C != 0 & data$Consumers_C == 0) {
            state = "Primary_Producers"
        } else if (data$Plants_C != 0 & data$Herbivores_C == 0 & data$Decomposers_C != 0 & data$Consumers_C != 0) {
            state = "H_extinct"
        } else if (data$Plants_C != 0 & data$Herbivores_C != 0 & data$Decomposers_C != 0 & data$Consumers_C == 0) {
            state = "C_extinct"
        } else {
            state = "no_organisms"
        }
    } else {
        if (data$Plants_C == 0 & data$Decomposers_C != 0) {
            state = "Terrestrial_extinct"
        } else if (data$Plants_C != 0 & data$Decomposers_C == 0) {
            state = "Aquatic_extinct"
        } else if (data$Plants_C != 0 & data$Herbivores_C != 0 & data$Decomposers_C != 0 & data$Consumers_C != 0 & data$Top_cons_H_C != 0 & data$Top_cons_C_C != 0) {
            state = "Coexistence"
        } else if (data$Plants_C != 0 & data$Herbivores_C == 0 & data$Decomposers_C != 0 & data$Consumers_C == 0 & data$Top_cons_H_C == 0 & data$Top_cons_C_C == 0) {
            state = "Primary_Producers"
        } else if (data$Plants_C != 0 & data$Herbivores_C != 0 & data$Decomposers_C != 0 & data$Consumers_C != 0 & data$Top_cons_H_C == 0 & data$Top_cons_C_C != 0) {
            state = "Top H extinct"
        } else if (data$Plants_C != 0 & data$Herbivores_C != 0 & data$Decomposers_C != 0 & data$Consumers_C != 0 & data$Top_cons_H_C != 0 & data$Top_cons_C_C == 0) {
            state = "Top C extinct"
        } else {
            state = "no_organisms"
        }
    }
    return(state)
}

Primary_production = function(state, param, consumers = F, colim = F, DC = F) {
    if (class(state)[1] == "tbl_df") {
        state = state[nrow(state), -1]
    }
    if (ncol(state) < 30) { # i.e. there is not the parameters
        state = cbind(state, matrix(unlist(param), ncol = length(param), nrow = 1))
        col_begin = ifelse(consumers, 19, 15)
        colnames(state)[col_begin:ncol(state)] = names(param)
    }

    limit = Get_limitation(state, param)

    if (limit$Limitation == "C-limited") {
        prod_aqua = state$eB * state$aBD * state$Detritus_A_C
    }
    if (limit$Limitation == "N-limited") {
        prod_aqua = state$aBN * state$Nitrogen_A_N
    }

    if (DC == T) {
        prod_terr = state$aP * state$Nitrogen_T_N
    } else {
        prod_terr = state$aP * state$Nitrogen_T_N * state$Plants_C
    }


    if (colim) {
        prod_aqua = state$eB * state$aBD * state$Detritus_A_C * state$aBN * state$Nitrogen_A_N
    }

    prod_tot_aq = state$aBD * state$Decomposers_C * state$Detritus_A_C + state$aBN * state$Decomposers_C * state$Nitrogen_A_N # for testing

    return(list(
        Terrestrial = prod_terr,
        Aquatic = prod_aqua, Aqua_tot = prod_tot_aq
    ))
}
