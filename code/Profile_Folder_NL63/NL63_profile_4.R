#### GENERAL ######
library(tidyverse)
library(readxl) ## read excel files

#### MODELS #####
library(deSolve) ## ODE solver
library(bbmle) ## likelihood analysis

#### Get job ID ####
ID <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()
print(ID)
#ID <- 1
#### Get rho ####
rho_95 <- 0.3912531

#### IVM With 5 groups ####
variation_infection_model_4 <- function(t, y, pars){ 
  # State variables
  P_0 <- y[1]
  
  M_1 <- y[2]
  P_1 <- y[3]
  
  M_2 <- y[4]
  P_2 <- y[5]
  
  M_3 <- y[6]
  P_3 <- y[7]
  
  M_4 <- y[8]
  P_4 <- y[9]
  
  
  # Parameter values
  rho_m <- pars["rho_m"] #maternal immunity waning stage
  w1 <- pars["w1"]
  w2 <- pars["w2"] 
  w3 <- pars["w3"] 
  rho_4 <- pars["rho_4"] 
  
  l43 <- pars["l43"]
  l32 <- pars["l32"]
  l21 <- pars["l21"]
  l10 <- pars["l10"]
  
  b0 <- pars["b0"]
  b1 <- pars["b1"]
  b2 <- pars["b2"]
  
  # calculated params 
  rho_3 <- w3*rho_4
  rho_2 <- w2*rho_3
  rho_1 <- w1*rho_2
  
  l42 <- b2*l43
  l41 <- b1*l42
  l40 <- b0*l41
  
  l31 <- b1*l32
  l30 <- b0*l31
  
  l20 <- b0*l21
  
  # rates in 
  seroconvert_into_1 <- l10*P_0
  seroconvert_into_2 <- l20*P_0 + l21*P_1
  seroconvert_into_3 <- l30*P_0 + l31*P_1 + l32*P_2
  seroconvert_into_4 <- l40*P_0 + l41*P_1 + l42*P_2 + l43*P_3
  
  maternal_revert_into_0 <- rho_m*(M_1 + M_2 + M_3 + M_4)
  
  # rates out 
  seroconvert_outof_0 <- l10*P_0 + l20*P_0 + l30*P_0 + l40*P_0
  seroconvert_outof_1 <- l21*P_1 + l31*P_1 + l41*P_1
  seroconvert_outof_2 <- l32*P_2 + l42*P_2
  seroconvert_outof_3 <- l43*P_3
  
  # Equations
  
  dP_0 <- -seroconvert_outof_0 + maternal_revert_into_0 + rho_1*P_1
  
  dM_1 <- -rho_m*M_1
  dP_1 <- +seroconvert_into_1  - seroconvert_outof_1 - rho_1*P_1 + rho_2*P_2
  
  dM_2 <- -rho_m*M_2
  dP_2 <- seroconvert_into_2 - seroconvert_outof_2 - rho_2*P_2 + rho_3*P_3
  
  dM_3 <- -rho_m*M_3
  dP_3 <- seroconvert_into_3 - seroconvert_outof_3 - rho_3*P_3 + rho_4*P_4
  
  dM_4 <- -rho_m*M_4
  dP_4 <- seroconvert_into_4 - rho_4*P_4
  
  
  # Return list of gradients 
  out <- c(dP_0, dM_1, dP_1, dM_2, dP_2, dM_3, dP_3, dM_4, dP_4)
  list(out) 
}


#### Log likelihood ####
loglik4_variation <- function(rho_m, 
                              rho,
                              l43, l32,
                              l21, l10, 
                              b0, b1,
                              b2,
                              sigma){
  
  paras <- c(rho_m = rho_m, w1 = 1,
             w2 = 1, w3 = 1,
             rho_4 = rho, 
             l43 = l43, l32 = l32,
             l21 = l21, l10 = l10, 
             b0 = b0, b1 = b1,
             b2 = b2,
             sigma = sigma)
  # print(pars)
  
  times <- seq(t_start,70, by = 0.125) #by = .125)
  
  serocatalytic_model_out <- ode(y = init, times = times, func = ODE_fun, parms = paras) %>% as.data.frame()
  
  ## double check the population math 
  df_out <- serocatalytic_model_out %>%
    pivot_longer(-time, names_to = "State", values_to = "Value") %>% 
    separate(State, into = c("Type", "Serostatus"), sep = "_") %>% 
    group_by(time) %>%
    mutate(pop = sum(Value)) %>%
    ungroup() %>% 
    group_by(Serostatus, Type, time, pop) %>% 
    summarise(.groups = "keep", Value = sum(Value)) %>% 
    mutate(Prop = Value/pop) %>%
    ungroup() %>% 
    rename(Group = Serostatus) %>%
    mutate(Group = as.numeric(Group)+1) %>%
    select(time, Prop, Type, Group) %>%
    inner_join(data %>% select(time = Age_con,Group, Type), by = c("time","Type", "Group")) %>% 
    arrange(time, Type, Group)
  
  ll <- -sum(dnorm(x=data$sero_prop,mean=df_out$Prop,sd=sigma,log=TRUE))
  ll = if_else(is.na(ll) | is.infinite(ll), 10^6, ll)
  
  return(ll)
}


#### Profile function ####
profile_variation_4 <- function(i){
  
  dt_profile <- read_csv("data_in/return_fits.csv") %>%
    filter(strain == "NL63_S1", rho == "rho_95", boundaries == "Bounded", waning == "Laddered", exposure_hypothesis == "Variation") %>%
    select(paras, Estimate, loglik) %>% 
    mutate(upper = unlist(unname(c(list(rho_m = 365,
                                        l43 = 12, l32 = 12, l21 = 12, l10 = 12,
                                        b0 = 8000, b1 = 8000, b2 = 8000, sigma = 5)))),
           lower = unlist(unname(c(list(rho_m = 0,
                                        l43 = 1/80, l32 = 1/80, l21 = 1/80, l10 = 1/80,
                                        b0 = 1, b1 = 1, b2 = 1, sigma = 0.001))))
    )
  
  
  data <- read_csv("data_in/data_all.csv") %>% filter(strain == "NL63_S1") %>% arrange(Age_con, Type, Group)
  
  ## initial conditions
  data_start <- data %>% filter(Age_con == 0)
  init <- c(P_0 = data_start$sero_prop[data_start$Group == 1], 
            M_1 = data_start$sero_prop[data_start$Group == 2], P_1 = 0,
            M_2 = data_start$sero_prop[data_start$Group == 3], P_2 = 0,
            M_3 = data_start$sero_prop[data_start$Group == 4], P_3 = 0,
            M_4 = data_start$sero_prop[data_start$Group == 5], P_4 = 0
  )
  
  ## globally assign init
  assign("init", init, envir = .GlobalEnv)
  assign("data", data, envir = .GlobalEnv)
  
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms.
  start2 = list(
    rho_m = dt_profile$Estimate[1], 
    l43 = dt_profile$Estimate[2], 
    l32 = dt_profile$Estimate[3], 
    l21 = dt_profile$Estimate[4], 
    l10 = dt_profile$Estimate[5], 
    b0 = dt_profile$Estimate[6], 
    b1 = dt_profile$Estimate[7], 
    b2 = dt_profile$Estimate[8], 
    sigma = dt_profile$Estimate[9]
  )
  upper_vec = list(rho_m = 365,  l43 = 12, l32 = 12, l21 = 12, l10 = 12,
                   b0 = 8000, b1 = 8000, b2 = 8000,  sigma = 5)
  lower_vec = list(rho_m = 0,  l43 = 1/80, l32 = 1/80, l21 = 1/80, l10 = 1/80,
                   b0 = 1, b1 = 1, b2 = 1,  sigma = 0.001)
  fixed_list2 = list(
    rho_m =  parm_test,
    l43 = parm_test,
    l32 = parm_test,
    l21 = parm_test,
    l10 = parm_test,
    b0 = parm_test,
    b1 = parm_test,
    b2 = parm_test,
    sigma = parm_test
  )
  
  startB <- start2[-i] ## remove the one we are fixing
  upper_listB <- upper_vec[-i]
  lower_listB <- lower_vec[-i]
  fixed_listB <- fixed_list2[i] ## make a list containing only the param we are fixing
  
  conf_threshold <- dt_profile$loglik[1] - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 1*max(dt_profile$Estimate[i],0.001)
  parm_test_right = dt_profile$Estimate[i] +  precision ## start at the true parameter value plus a step
  upper_limit_var <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    parm_test_right = dt_profile$Estimate[i]
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      if (parm_test_right >= upper_limit_var){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          rho_m =  parm_test_right,
          l43 = parm_test_right,
          l32 = parm_test_right,
          l21 = parm_test_right,
          l10 = parm_test_right,
          b0 = parm_test_right,
          b1 = parm_test_right,
          b2 = parm_test_right,
          sigma = parm_test_right
        )
        fixed_listB <- fixed_listB[i]
        
        
        ## run mle and grab the new log likelihood
        print(parm_test_right)
        m1_prof = mle2(minuslogl = loglik4_variation,
                       start = startB,
                       data = data,
                       fixed = c(fixed_listB, rho = rho_95),
                       method="L-BFGS-B",
                       upper= upper_listB,
                       lower = lower_listB,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        right_ll <- -m1_prof@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        
        new_ll_right <- right_ll
      }
      
      ## If it equals the confidence threshold, stop
      if (round(new_ll_right, 4) == round(conf_threshold, 4)){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  
  ## First, set the starting steps away from the parameter
  precision = 1*max(dt_profile$Estimate[i],0.001)
  parm_test_left = dt_profile$Estimate[i]  - precision ## start at the true parameter value minus first step
  lower_limit_var <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$lower[i],4)){
    parm_test_left = dt_profile$Estimate[i]
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit_var){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        
        
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          rho_m =  parm_test_left,
          l43 = parm_test_left,
          l32 = parm_test_left,
          l21 = parm_test_left,
          l10 = parm_test_left,
          b0 = parm_test_left,
          b1 = parm_test_left,
          b2 = parm_test_left,
          sigma = parm_test_left
        )
        fixed_listB <- fixed_listB[i]
        
        #print(parm_test_left)
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        print(parm_test_left)
        m1_prof = mle2(minuslogl = loglik4_variation,
                       start = startB,
                       data = data,
                       fixed = c(fixed_listB, rho = rho_95),
                       method="L-BFGS-B",
                       upper= upper_vec,
                       lower = lower_vec,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        left_ll <- -m1_prof@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        new_ll_left <- left_ll
        
      }
      
      ## If it equals the confidence threshold, stop
      if (round(new_ll_left,4) == round(conf_threshold,4)){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  out_data <- tibble(param = dt_profile$paras[i],
                     loglik = dt_profile$loglik[1],
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("data_out/", dt_profile$paras[i], ".csv"))
  
}


ODE_fun <- variation_infection_model_4
t_start = 0 # set start time

profile_variation_4(ID)




