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

#### Binary model ####
serocatalytic_2_model <- function(t, y, pars){ 
  # State variables
  P_0 <- y[1]
  M_1 <- y[2]
  P_1 <- y[3]
  
  # Fit parameter values
  rho_m <- pars["rho_m"] #maternal immunuty wanning stage
  rho_1 <- pars["rho_1"] # maximum waning of immmunity
  lambda_0 <- pars["lambda_0"]
  
  # Equations
  dP_0 <- rho_1*P_1 - lambda_0*P_0 + (rho_m*M_1)
  
  dM_1 <- -rho_m*M_1
  dP_1 <- + lambda_0*P_0  - rho_1*P_1
  
  # Return list of gradients 
  out <- c(dP_0, dM_1, dP_1)
  list(out) 
}


#### Log likelihood ####
loglik2 <- function(rho_m,
                    rho,
                    lambda_0, sigma){
  
  paras <- c(rho_m = rho_m,
             rho_1 = rho, 
             lambda_0 = lambda_0, 
             sigma = sigma)
  #print(paras)
  
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
    mutate(Group = as.numeric(Group) + 1) %>%
    select(time, Prop, Type, Group) %>%
    inner_join(data %>% 
                 select(time = Age_con, Group, Type, sero_prop), 
               by = c("time", "Type", "Group")) %>% 
    arrange(time, Type, Group) 
  
  ll <- -sum(dnorm(x=df_out$sero_prop,mean=df_out$Prop,sd=sigma,log=TRUE))
  ll = if_else(is.na(ll) | is.infinite(ll), 10^6, ll)
  
  return(ll)
}


#### Profile function ####
profile_2 <- function(i){
  
  dt_profile <- read_csv("data_in/return_fits.csv") %>%
    filter(strain == "a_229E_S1", rho == "rho_95", boundaries == "Bounded", waning == "Direct", exposure_hypothesis == "Binary") %>%
    select(paras, Estimate, loglik) %>% 
    mutate(upper = unlist(unname(c(list(rho_m = 365, lambda_0 = 12, sigma = 5)))),
           lower = unlist(unname(c(list(rho_m = 0, lambda_0 = 1/80, sigma = 0.001))))
    )
  
  
  data <- read_csv("data_in/data_all2.csv") %>% filter(strain == "a_229E_S1") %>% arrange(Age_con, Type, Group)
  
  ## initial conditions
  data_start <- data %>% filter(Age_con == 0)
  init <- c(P_0 = data_start$sero_prop[data_start$Group == 1], 
            M_1 = data_start$sero_prop[data_start$Group == 2], P_1 = 0
  )
  
  ## globally assign init
  assign("init", init, envir = .GlobalEnv)
  assign("data", data, envir = .GlobalEnv)
  
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms.
  start2 = list(
    rho_m = dt_profile$Estimate[1], 
    lambda_0 = dt_profile$Estimate[2], 
    sigma = dt_profile$Estimate[3]
  )
  upper_vec = list(rho_m = 365, lambda_0 = 12, sigma = 5)
  lower_vec = list(rho_m = 0, lambda_0 = 1/80, sigma = 0.001)
  fixed_list2 = list(
    rho_m =  parm_test,
    lambda_0 = parm_test,
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
          lambda_0 = parm_test_right,
          sigma = parm_test_right
        )
        fixed_listB <- fixed_listB[i]
        
        
        ## run mle and grab the new log likelihood
        print(parm_test_right)
        m1_prof = mle2(minuslogl = loglik2,
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
          lambda_0 = parm_test_left,
          sigma = parm_test_left
        )
        fixed_listB <- fixed_listB[i]
        
        #print(parm_test_left)
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        print(parm_test_left)
        m1_prof = mle2(minuslogl = loglik2,
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


ODE_fun <- serocatalytic_2_model
t_start = 0 # set start time

profile_2(ID)




