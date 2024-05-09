#libraries 
library(dplyr)
library(tidyverse)
library(BOIN)
library(mvtnorm)
library(coda)
library(rriskDistributions)
library(nimble)

#---------------------------------------------------------------------------------------------------
# weibull_parms()

weibull_parms <- function(pDLT, DLT_time){
  #x1 = 0.5*(DLT_time) x2 = DLT_time
  shape <- (log(-log(1 - pDLT)) - log(-log(1 - 0.5*pDLT)))/(log(DLT_time) - log(0.5*DLT_time) )
  scale <- (0.5*DLT_time)/((-log(1 - 0.5*pDLT))^(1/shape))
  
  return(c(shape, scale))
}

#cohort_patient()
cohort_patient <- function(cohort, doses_info, current_dose, run, pts, previous_time, last_arrival_time = 0, last_lag = 0, i_simulation = 0){
  
  
  seed <- 1234 + run + i_simulation
  set.seed(seed)
  
  #take the parms for the weibull DLT according to the current_dose
  parms_DLT <- doses_info %>% filter(Dose == current_dose) %>% select(shape = Shape_DLT, scale = Scale_DLT)
  n_call <- 0
  for (i in 1:cohortsize){
    n_call <- n_call + 1
    if (pts == 0){
      pts <- pts + 1
      time_DLT <- round(rweibull(1, shape = parms_DLT$shape, scale = parms_DLT$scale), 4)
      cohort[pts, ] <- c(run, pts, 'C', current_dose, 0, time_DLT, ifelse(time_DLT <= DLT_time, 1, 0), 0, min(time_DLT, DLT_time), time_DLT)
      
    } else {
      pts <- pts + 1
      lag_arrival <- round(ifelse((last_arrival_time != 0 && n_call == 1), last_lag ,rexp(1, rate = lambda)), 4)
      time_arrival <- round(ifelse((last_arrival_time != 0 && n_call == 1), last_arrival_time, previous_time + lag_arrival ), 4)
      time_DLT <- round(rweibull(1, shape = parms_DLT$shape, scale = parms_DLT$scale), 4)
      cohort[pts, ] <- c(run, pts, 'C', current_dose, time_arrival, time_DLT, ifelse(time_DLT <= DLT_time, 1, 0), lag_arrival, min((DLT_time + time_arrival), (time_DLT + time_arrival)), (time_DLT + time_arrival) )
      previous_time <- time_arrival 
    }
    
    #the limit time is the previous time from which the next patients can be introduced 
    limit_time <- cohort %>% filter(Group == 'C')  %>% filter(Run == run) %>% summarise(limit_time = max(Limit_time))
    limit_time <- as.numeric(limit_time$limit_time)
    
  }
  
  
  return(list('cohort' = cohort, 'pts' = pts, 'previous_time' = previous_time, 'limit_time' = limit_time))
  
}

#backfill_patients()
backfill_patients <- function(cohort, run, pts, previous_time, limit_time, n_max, i_simulation = 0){
    
    seed <- 555 + run + i_simulation
    set.seed(seed)
    #generate potentially untill full of pts 
    while(pts < n_max){
      lag_arrival <- round(rexp(1, rate = lambda), 4)
      time_arrival_backfill <- round(previous_time + lag_arrival, 4)
      #if the pts arrives after limit_time (point in time-line for cohort)
      if(time_arrival_backfill >= limit_time){ break }
      pts <- pts + 1
      cohort[pts, ] <- c(run, pts, 'B', NA, time_arrival_backfill, NA, NA, lag_arrival, NA, NA)
      previous_time <- time_arrival_backfill 
      #previous_time keeps track of the last point in the time line for the backfill --> can be not returned. Not necessary for next cohort and next backfill
        
    }

  
  return(list('cohort' = cohort, 'pts' = pts, 'previous_time' = previous_time, 'last_arrival' = time_arrival_backfill, 'last_lag' = lag_arrival))
}

#chack_individual_posterior

check_individual_posterior <- function(pts, fail, dose, info_probs){
  
  alpha <- info_probs %>% filter(Dose == dose) %>% select(alpha)
  alpha <- as.numeric(alpha$alpha) + fail
  
  beta <- info_probs %>% filter(Dose == dose) %>% select(beta)
  beta <- as.numeric(beta$beta) + pts - fail 
  
  EWOC <- ifelse(pbeta(0.35, shape1 = alpha, shape2 = beta , lower.tail =  F) >= 0.25, F, T)
  
  return(EWOC)
}


#open_close_EWOC

open_close_EWOC <- function(cohort, doses_info, current_dose, time_arrival = 1000, run, target = 0.3, i_simulation = 0){
  
  info_probs <- data.frame('Dose' = doses_info$Dose, 'alpha' = rep(NA, nrow(doses_info), 'beta' = rep(NA, nrow(doses_info))))
  
  prob <- data.frame(matrix(ncol = length(doses_info$Dose), nrow = 9000)) #nrow = iterations - burnin
  
  #directly call the MCMC
  results <- MCMC_adaptive_EWOC(cohort, doses_info, time_arrival, i_simulation = i_simulation, run = run) #iterations = 10000, burnin = 1000
  #results_nimble <- MCMC_nimble(cohort, doses_info, time_arrival, i_simulation = i_simulation, run = run)
  
  #compare with Nimble results 
  #plot(as.mcmc(results$betas[, 1]))
  #plot(as.mcmc(results$betas[, 2]))
  #plot(as.mcmc(results_nimble$betas))
  
  #compute the probabilities (exp the beta_1 as you have samples of log(beta_1))
  for (i in 1:length(doses_info$Dose) ){
    logit_prob <- results$betas[, 1] + results$betas[, 2]*log(doses_info$Dose[i]/reference_dose)
    prob[, i] <- exp(logit_prob)/(1 + exp(logit_prob))
    #hist(prob[, i])
  }
 
  #compute the alpha and beta parsm from the quantiles of the posteriors and store them to be used later 
  info_probs <- compute_parameters(info_probs, prob)
  names(info_probs) <- c('Dose', 'alpha', 'beta')
  
  #compute the dose index that maximizes the posterior for the targeted interval while controlling for the excessive toxicity interval   
  index_current_dose <- which(doses_info == current_dose)
  max_prob <- evaluate_posterior(prob, limit_dose = index_current_dose)
  #print(c('index', max_prob))
  choice <- ifelse(is.na(max_prob), NA, doses_info$Dose[max_prob])
 
  return(list('choice' = choice, 'diagnostics' = results$diagnostics, 'info_probs' = info_probs))
  
}

#maximum_open_EWOC()

maximum_open_EWOC <- function(cohort, time_pts_backfill, current_dose, doses_info, info_probs, run, last_time){
  
  current_backfill_dose <- doses_info %>% filter(State == 1)  %>% filter(Dose < current_dose) %>% arrange(Dose) %>% tail(1)
  current_backfill_dose <- current_backfill_dose$Dose
        
  #valid dose?
  if (identical(current_backfill_dose, numeric(0))) {
    #print('No length')
    current_backfill_dose <- NA
    return(list('backfill_dose' = current_backfill_dose))
  }
  
  #is this still okay? check the dose-individual posteriors for the current backfilldose and for the cohort 
  cohort_data <- cohort %>% filter(Dose == current_dose) %>% filter(Limit_time <= time_pts_backfill) %>% filter(Limit_time >= last_time) %>% filter(Run == run) %>% summarise('DLT' = sum(DLT), 'Pts' = n())
  backfill_data <- cohort %>% filter(Dose == current_backfill_dose) %>% filter(Limit_time <= time_pts_backfill) %>% filter(Limit_time >= last_time) %>% filter(Run == run) %>% summarise('DLT' = sum(DLT), 'Pts' = n())
  
  cohort_check <- check_individual_posterior(fail = cohort_data$DLT, pts = cohort_data$Pts, dose = current_dose, info_probs = info_probs)
  backfill_check <- check_individual_posterior(fail = backfill_data$DLT, pts = backfill_data$Pts, dose = current_backfill_dose, info_probs = info_probs)
  
  #if both are okay then the current backfilldose is the one proposed. Otherwise we run the MCMC to find the next optimal dose that is lower than the current dose. If none can be find then the patient is not assigned to a backfill dose at all 
  
  if(isFALSE(cohort_check) || isFALSE(backfill_check)){
    
   results <- open_close_EWOC(cohort, doses_info, current_dose = current_dose, time_arrival = time_pts_backfill, run = run, i_simulation = 0)

   #update the doses_info so at the next stage we have proper info
   
   for (i in 1: length(doses_info$Dose)) {
     #print(c('dose state', doses_info$Dose[i],doses_info$State[i]))
     
     if(is.na(results$choice)){
       #if no choice then all open ones set to 0
       doses_info$State[i] <- ifelse(doses_info$State[i] == 1, 0, doses_info$State[i])
     } else {
       #if ok choice then decide according to the dose level
       if (doses_info$State[i] == 1 || doses_info$State[i] == 0 ){
         if (doses_info$Dose[i] <= results$choice){
           doses_info$State[i] <- 1 
         } else {
         doses_info$State[i] <- 0}
       }
      } 
    }

   #imagine no results 
   
   if(is.na(results$choice)){
     current_backfill_dose <- NA
     return(list('backfill_dose' = current_backfill_dose, 'info_probs' = results$info_probs, 'diagnostics' = results$diagnostics, 'doses_info' = doses_info))
   } 
   
   #ensure that the backfill dose is always lower than the current cohort dose
   #backup_doses <- doses_info %>% filter(Dose <= current_dose) %>% select(doses = Dose)
   #choice <- ifelse(results$choice >= current_dose, backup_doses$doses[(length(backup_doses$doses) - 1)], results$choice )
   #print(choice)
   
   #imagine the selected dose (after the change. This after the change as we exclude some, while before we were closing or opening based on the true choice) can not be used since too many pts or not response, then do no put the pts
   
  index_dose <- doses_info %>% filter(Dose == results$choice) %>% select(State)
   #print(c('index dose', index_dose$State))
   if (index_dose$State == -1 || index_dose$State == 2){
     current_backfill_dose <- NA
   } else {current_backfill_dose <- results$choice}
   
   
    return(list('backfill_dose' = current_backfill_dose, 'info_probs' = results$info_probs, 'diagnostics' = results$diagnostics, 'doses_info' = doses_info))
    
   }
  
  return(list('backfill_dose' = current_backfill_dose))
        
}

#logposterior()

logposterior <- function(data, betas){
  
  #define the prior parameters as in the BRLM+BF paper (6 doses)
  prior_mean <- c(log(0.5), 0)
  prior_var <- matrix(c(4, 0 ,0, 1), nrow = 2)
  
  #vector of probabilities based on the betas and on the doses with the model dependign of beta0 and beta1 (transform with exp)
  logit_p <- betas[1] + exp(betas[2])*log(data$Dose/reference_dose)
  p <- 1/(1 + exp(-logit_p))

  likelihood <- sum(dbinom(data$DLT, data$Pts, p, log = T))
  prior <- dmvnorm(betas, mean = prior_mean, sigma = prior_var, log = T)
  
  #the det(jacobian) is the exp(beta2) that multiplies the posterior. Since work in log then it is the beta2  
 
  return(sum(sum(likelihood, prior), betas[2]) )
  
}

#MCMC
MCMC <- function(cohort, doses_info, time_arrival, i_simulation = 0, run, iterations = 10000, burnin = 1000){
  
  #set the seed according to the run and to the simulation number. Ok global variables  
  set.seed((1234 + run + i_simulation))
  
  #define the sd factor for the acceptance rate
  sd <- 2.4
  
  #take the values in a summary format
  dlt <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(freq = sum(DLT))
  pts <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(pts = n())
  doses <- dlt$Dose #select the appropriate doses
  
  data <- data.frame('DLT' = dlt$freq, 'Pts' = pts$pts, 'Doses' = dlt$Dose )
  
 
  beta0_start <- 1
  beta1_start <- 0
  sigma_start <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

  #implement MCMC calling the logpost function 
  
  #parameters 
  n_accept <- 0
  beta_parms <- matrix(nrow = iterations, ncol = 2)
  beta_parms[1, ] <- c(beta0_start, beta1_start)
  
 for( i in 2:iterations){
  
    beta_proposal <- rmvnorm(n = 1, mean = c(1, 0), sigma = (sd^2)*sigma_start) #sampling of beta0 and log(beta1)
     
    log_now <- logposterior(data, beta_proposal)
    log_previous <- logposterior(data, beta_parms[(i - 1), ])
    log_ratio <- log_now - log_previous
    
    if(is.na(log_ratio)){
      print('change')
      beta_parms[i, ] <- beta_parms[(i-1), ]
      #go to the next iteration
      next
    }
  
    u <- runif(1, min = 0, max = 1)
    if(log_ratio >= 0 || u <= exp(log_ratio)){
      beta_parms[i, ] <- beta_proposal
      n_accept <- n_accept + 1
    } else {
      beta_parms[i, ] <- beta_parms[(i-1), ]
    }
  }
  
  #trim the dataset for the burnin 
  beta_parms <- beta_parms[-c(1:burnin), ]
  
  #run and collect the diagnostics: Geweke and acceptance rate 
  diagnostics <- c(round(n_accept/(iterations - 1), digits = 2) * 100, geweke.diag(as.mcmc(beta_parms[, 1]))$z,  geweke.diag(as.mcmc(beta_parms[, 2]))$z )

  beta_parms[, 2] <- exp(beta_parms[, 2])
  #need for the parms and for the diagnostics
  return(list('betas' = beta_parms, 'diagnostics' = diagnostics))
  
}

#MCMC_adaptive_EWOC

MCMC_adaptive_EWOC <- function(cohort, doses_info, time_arrival, i_simulation = 0, run, iterations = 10000, burnin = 1000, delta = 0.01){
  
  #set the seed according to the run and to the simulation number. Ok global variables  
  set.seed((1234 + run + i_simulation))
  

  #take the values in a summary format
  dlt <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(freq = sum(DLT))
  pts <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(pts = n())
  doses <- dlt$Dose #select the appropriate doses
  
  data <- data.frame('DLT' = dlt$freq, 'Pts' = pts$pts, 'Doses' = dlt$Dose )
  print(data)

  mean_current <- c(1, 1)
  beta0_start <- 1 
  beta1_start <- 1
  sigma_current <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  
  #implement MCMC calling the logpost function 
  
  #parameters 
  n_accept <- 0
  beta_parms <- matrix(nrow = iterations, ncol = 2)
  colnames(beta_parms) <- c('beta0', 'beta1')
  beta_parms[1, ] <- c(beta0_start, beta1_start)
  
  #define the sd factor for the acceptance rate according to Andrieu and Thoms 
  lambda_proposed <- log((2.38^2)/2)
  
  for( i in 0:(iterations-2) ){
    
    beta_proposal <- rmvnorm(n = 1, mean = mean_current, sigma = exp(lambda_proposed)*sigma_current) #sampling of beta0 and log(beta1)
     
    log_now <- logposterior(data, beta_proposal)
    log_previous <- logposterior(data, beta_parms[(i - 1 + 2), ])
    log_ratio <- log_now - log_previous
   
    if(is.na(log_ratio)){
      
      print('change')
      #print('notaccepted')
      beta_parms[(i + 2), ] <- beta_parms[(i-1 + 2), ]
      
      #go to the next iteration
      next
      
    }
    
    u <- runif(1, min = 0, max = 1)
    if(log_ratio >= 0 || u <= exp(log_ratio)){
      beta_parms[(i + 2), ] <- beta_proposal
      n_accept <- n_accept + 1
    } else {
      beta_parms[(i + 2), ] <- beta_parms[(i - 1 + 2), ]
    }
    
    #updating 
    
    #weights
    gamma <- 1/(1 + (i))^delta  #1/(1 + (i - 1))
    
    #update other parms 
    #lambda tuning on the logscale
    lambda_proposed <- gamma*(min(exp(log_ratio), 1) - 0.24) + lambda_proposed
    #sigma varcov
    sigma_current <- sigma_current + gamma*((beta_parms[(i + 2), ] - mean_current)%*%t(beta_parms[(i + 2), ] - mean_current) - sigma_current ) + 0.001 #can i add a small quantity here?
    
    #mean vector
    mean_current <- mean_current + gamma*(beta_parms[(i + 2), ] - mean_current)    
  }
  
  #trim the dataset for the burnin 
  beta_parms <- beta_parms[-c(1:burnin), ]
  
  #run and collect the diagnostics: Geweke and acceptance rate 
  diagnostics <- c(round(n_accept/(iterations - 1), digits = 2) * 100, geweke.diag(as.mcmc(beta_parms[, 1]))$z,  geweke.diag(as.mcmc(beta_parms[, 2]))$z )

 
  #need for the parms and for the diagnostics
  beta_parms[, 2] <- exp(beta_parms[, 2])
  return(list('betas' = beta_parms, 'diagnostics' = diagnostics))
  
}

#MCMC_nimble

MCMC_nimble <- function(cohort, doses_info, time_arrival, i_simulation = 0, run, iterations = 10000, burnin = 1000, delta = 0.01){
  
  #set the seed according to the run and to the simulation number. Ok global variables  
  set.seed((1234 + run + i_simulation))
  
  data <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) 
  
  
  model1 <- nimbleCode(
  {for (i in 1:pts){
    
    logit(pi[i]) <- beta[1] + exp(beta[2])*Dose[i]
    et[i] ~ dbern(pi[i])
    
  }
  
  beta[1] ~ dnorm(log(0.5), sd = 2)
  beta[2] ~ dnorm(0, sd = 1)
} )
  
  
  N <- nrow(data)
  my.data <- list(et = data$DLT, Dose = data$Dose)
  my.constants <- list(pts = N)

  my.inits <- list(list(beta = c(0,0)))

  parameters <- c('beta')

  model.sim <- nimbleMCMC(code = model1, 
                        data = my.data, 
                        constants = my.constants, 
                        inits = my.inits, 
                        monitors = parameters, 
                        niter = 10000,
                        nburnin = 1000,
                        nchains = 1, 
                        thin = 1, 
                        samplesAsCodaMCMC = T)
  
  return(list('betas' = model.sim, 'diagnostics' = c(NA, NA, NA)))
}

#evaluate_posterior()
evaluate_posterior <- function(prob, limit_dose = NA){
  
  targeted_interval <- array()
  excessive_interval <- c()
  excessive <- T
  
  for(i in 1:ncol(prob)){
    
    values <- prob[, i]
    targeted_interval[i] <- mean(values < 0.35) - mean(values <= 0.20)
    excessive_interval[i] <- mean(values >= 0.35)
 
  }
  
  #print(c('targeted', round(targeted_interval, 2)))
  #print(c('excessive', round(excessive_interval, 2)))
  
  targeted_interval <-  round(targeted_interval, 2)
  excessive_interval <-  round(excessive_interval, 2)
  #see open_close_EWOC. to select doses up to the current cohort dose not included 
  
  targeted_interval <- ifelse(is.na(targeted_interval), 0, targeted_interval)
  
  if (!is.na(limit_dose)){targeted_interval <- targeted_interval[1:(limit_dose - 1)]
  }
  
  while(isTRUE(excessive)){
    
    #imagine that you have already looped over all the vector and no dose is fine so then all the values are negative: no dose is okay so NA and exit
    if (sum(targeted_interval > 0) == 0){
      #print(c('ended ti', targeted_interval))
      dose_index <- NA
      break
    }
    
    dose_index <- which(targeted_interval == max(targeted_interval))
    dose_index <- dose_index[length(dose_index)]
    excessive <- ifelse(excessive_interval[dose_index] >= 0.25, T, F) #is the same dose ok for the toxicity? If yes, then exit, otherwise loop again to find the next maximizing 
    targeted_interval[dose_index] <- -targeted_interval[dose_index] #exclude from the next possible iteration
    
  }
  
  #return the dose index that has been found 
  return(dose_index)
}

#compute_parameters()

compute_parameters <- function(info_probs, probs){
  
  for(i in 1:ncol(probs)){
    
    qval <- as.numeric(quantile(probs[, i], c(0.025, 0.5, 0.975)))
    print(qval)
    parms <- get.beta.par(q = qval, show.output = F, plot = F)
    alpha <- as.numeric(parms[1])
    beta <- as.numeric(parms[2])
    
    info_probs[i, 2:3] <- c(round(alpha, 2), round(beta, 2))

  }
  
 return(info_probs)
  
}

#decision_EWOC()
decision_EWOC <- function(cohort, doses_info, time_arrival = 1000, run, target = 0.3, i_simulation = 0){
  
  info_probs <- data.frame('Dose' = doses_info$Dose, 'alpha' = rep(NA, nrow(doses_info), 'beta' = rep(NA, nrow(doses_info))))
  prob <- data.frame(matrix(ncol = length(doses_info$Dose), nrow = 9000)) #nrow = iterations - burnin
  
  #directly call the MCMC
  results <- MCMC_adaptive_EWOC(cohort, doses_info, time_arrival, i_simulation = i_simulation, run = run) #iterations = 10000, burnin = 1000

  #compare the results with Nimble 
  #plot(as.mcmc(results$betas[, 1]))
  #plot(as.mcmc(results$betas[, 2]))
  
  results_nimble <- MCMC_nimble(cohort, doses_info, time_arrival, i_simulation = i_simulation, run = run)
  #plot(as.mcmc(results_nimble$betas))
  
  #compute the probabilities (exp the beta_1 as you have samples of log(beta_1) --> changed to have directly beta1 look at MCMC)
  for (i in 1:length(doses_info$Dose) ){
    logit_prob <- results$betas[, 1] + results$betas[, 2]*log(doses_info$Dose[i]/reference_dose)
    prob[, i] <- exp(logit_prob)/(1 + exp(logit_prob))
    #hist(prob[, i])
  }
 
  
  
  #compute the alpha and beta parsm from the quantiles of the posteriors and store them to be used later 
  info_probs <- compute_parameters(info_probs, prob)
  names(info_probs) <- c('Dose', 'alpha', 'beta')
  
  #compute the dose index that maximizes the posterior for the targeted interval while controlling for the excessive toxicity interval   
  max_prob <- evaluate_posterior(prob)
  choice <- ifelse(is.na(max_prob), NA, doses_info$Dose[max_prob])
  print(c('choice', choice))
  
  return(list('choice' = choice, 'diagnostics' = results$diagnostics, 'info_probs' = info_probs))
  
}

#check_n_cap()
check_n_cap <- function(cohort, doses_info){
  
  state_update <- doses_info$State
  doses <- doses_info$Dose
  close_doses <- cohort %>% group_by(Dose) %>% summarise(tot_dose = n() ) %>% filter(tot_dose >= n_cap) %>% select(Dose)
  for(dose in close_doses$Dose){
    state_update[which(doses == dose)] <- 2
  }
  
  doses_info$State <- state_update
  
  return(doses_info)
  
}

#transform_data()
transform_data <- function(cohort){
  
  suppressWarnings(cohort2 <- transform(cohort, Run = as.numeric(Run), Pts = as.numeric(Pts) , Dose = as.numeric(Dose), Time_arrival = as.numeric(Time_arrival), Time_DLT = as.numeric(Time_DLT), DLT = as.numeric(DLT), Lag = as.numeric(Lag), Limit_time = as.numeric(Limit_time), Sum_times = as.numeric(Sum_times) ))
  
  return(cohort2)
  
}

#transform_data_doses()
transform_data_doses <- function(doses_info){
  
  suppressWarnings(cohort2 <- transform(doses_info, Dose = as.numeric(Dose), Shape_DLT = as.numeric(Shape_DLT) , Scale_DLT = as.numeric(Scale_DLT), Prob_resp = as.numeric(Prob_resp), Shape_resp = as.numeric(Shape_resp), Scale_resp = as.numeric(Scale_resp), State = as.numeric(State)))
  
  return(cohort2)
  
}


#compute_BLRMBF_EWOC()
compute_BFBLRM_EWOC <- function(cohort, doses_info, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda, target = 0.3 ,i_simulation = 0){
  
  run <- 0 
  pts <- 0
  n_early_stop <- 0
  n_stop_reached <- 0
  last_closed_dose <- NA
  limit_time <- 0
  current_dose <- min(doses) #start from the lowest dose
  last_arrival <- 0
  last_lag <- 0
  doses <- doses_info$Dose
  diagnostics <- matrix(ncol = 3)
    
  #initiate the while loop that should be stopped asap we reach: 1) n_max, 2) overtoxicity at minimal dose, 3) n_stop + stay decision 

  while (pts <= n_max){
    
    run <- run + 1
    
    #no backfilling
    if (sum(is.na(doses_info$State)) != 0 ){
        
      update <- cohort_patient(cohort, doses_info, current_dose, run, pts, previous_time = limit_time)
      cohort <- transform_data(update$cohort)
      pts <- update$pts
      previous_time <- update$previous_time
      limit_time <- update$limit_time
      #check the number of pts
      if(pts >= (n_max - cohortsize + 1)){break}
      
      #compute the decision and save the diagnostics 
      results <- decision_EWOC(cohort, doses_info, target = target, run = run, i_simulation = i_simulation)
      diagnostics <- rbind(diagnostics, results$diagnostics)
      info_probs <- results$info_probs
      #deciison + safety rules at section 3.3
      next_dose <- results$choice
      
      #if no dose is okay then we have to stop (so either it is the minimum or not)
      if (is.na(next_dose) && pts > 3){
        n_early_stop <- n_early_stop + 1
        break
      }
       
      #if not break then continue
       if(is.na(next_dose)){
        next_dose <- current_dose
      }
      
      #check for SUfficient information (same as 3cohorts + 'STAY')
      if (next_dose == current_dose){
        Pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
        if(Pts$pts >= n_stop){
          n_stop_reached <- n_stop_reached + 1
          break
        }
      }
      
      #if the next dose for the cohort is not the minimum dose then the backfill can start; if the current dose is the maximum one then the backfill option is not explored for the current dose level 
      backfill <- ifelse(next_dose != min(doses), ifelse(current_dose != max(doses), T, F), F)
      
            
      #possibility of backfilling
      if(backfill == T){
        parms_resp <- doses_info %>% filter(Dose == current_dose) %>% select(shape = Shape_resp, scale = Scale_resp)
        n_resp <- 0
      for (i in (pts-cohortsize):(pts - 1) ){
        resp <- rweibull(1, shape = parms_resp$shape, scale = parms_resp$scale)
        resp_i <- ifelse(resp < min((DLT_time + as.numeric(cohort[(i+1), ]$Time_arrival)), limit_time), 1, 0)
        n_resp <- n_resp + resp_i
      }
      
      #index <- which(doses == current_dose)
      if (n_resp >= 1){
        state <- ifelse(doses_info$Dose >= current_dose, 1, ifelse(doses_info$Dose < current_dose, -1, 0))
        doses_info$State <- state
        } 
      } 
  } else {
    
    #suppose the backfill is open so then we can work with this 
    print(doses_info)
  
    #vector of backfill has been explored already (no NA)
    
    # cohort pts --> cohortsize added
    update <- cohort_patient(cohort, doses_info, current_dose, run, pts, previous_time = limit_time, last_arrival_time = last_arrival, last_lag = last_lag)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    #previous_time for backfill adn limit_time for cohort
    previous_time <- update$previous_time
    limit_time <- update$limit_time
    #check the number of pts 
    if(pts >= n_max){break}


    #simulate arrival times 
    update <- backfill_patients(cohort, run, pts, previous_time = previous_time, limit_time = limit_time, n_max)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    #save last arrivals
    last_arrival <- update$last_arrival
    last_lag <- update$last_lag
    
    
    #decide the dose level for each of the arrived patients 
    temp_pts <- cohort %>% filter(Run == run) %>% filter(Group == 'B') %>% select(pts = Pts)
    temp_pts <- temp_pts$pts
    
    
    #how to check if the backfill dose is still ok? use MCMC in case you have not good a posteriori cohort and backfill (check with individual doses posteriori)
    discarded <- 0
    last_time <- cohort %>% filter(Group == 'C') %>% filter(Run == run) %>% head(1) %>%  select(time = Time_arrival)
    for (back_pts in temp_pts){
      #actual time of arrival to be considered to check if the current dose for backfill is stll okay 
      time_pts_backfill <- cohort %>% filter(Pts == back_pts) %>% select(time = Time_arrival)
      update_backfill <- maximum_open_EWOC(cohort, time_pts_backfill$time, current_dose, doses_info, info_probs, run = run, last_time = last_time$time)
      
      current_maximum_backfill_dose <- update_backfill$backfill_dose
      print(current_maximum_backfill_dose)
      
      print(update_backfill$doses_info)
      
      if(!is.null(update_backfill$doses_info)){
        doses_info <- transform_data_doses(update_backfill$doses_info)
      }
      
      print(update_backfill$diagnostics)
      
      if(!is.null(update_backfill$diagnostics)){
        diagnostics <- rbind(diagnostics, update_backfill$diagnostics)
      }
      
      print(update_backfill$info_probs)
      
      if(!is.null(update_backfill$info_probs)){
        info_probs <- update_backfill$info_probs 
      }
  

      #no backfill doses is open by arrival of the pts? then we discard the patient 
      if (is.na(current_maximum_backfill_dose)){
        discarded <- discarded + 1 
        next
      }
      
      #assign the dose and compute the time_DLT
       parms_backfill_DLT <- doses_info %>% filter(Dose == current_maximum_backfill_dose) %>% select(shape = Shape_DLT, scale = Scale_DLT)
      time_DLT <- round(rweibull(1, shape = parms_backfill_DLT$shape, scale = parms_backfill_DLT$scale), 4)
      cohort[back_pts, ]$Dose <- current_maximum_backfill_dose
      cohort[back_pts, ]$Time_DLT <- time_DLT
      cohort[back_pts, ]$DLT <- ifelse(time_DLT <= DLT_time, 1, 0)
      cohort[back_pts, ]$Limit_time <- min((DLT_time + cohort[pts,]$Time_arrival), (time_DLT + cohort[pts, ]$Time_arrival) )
      cohort[back_pts, ]$Sum_times <- time_DLT+ cohort[pts, ]$Time_arrival
      last_time <- cohort %>% filter(Group == 'B') %>% filter(!is.na(Dose)) %>% filter(Run == run) %>% tail(1) %>%  select(time = Time_arrival)
      
      #after adding the new pts, check if the states are all ok in terms of n_cap
      doses_info <- check_n_cap(cohort, doses_info)
    }
   
    #remove all the non-assigned patients entry and re-evaluate their number
    cohort <- cohort %>% filter(!is.na(Dose)) 
    if ((discarded - length(temp_pts)) != 0){
      count <- cohort %>% filter(Group == 'C') %>% tail(1) %>% select(Pts)
      count <- count$Pts
      for (value in (discarded - length(temp_pts) + 1):0){
        count <- count + 1
        cohort[(nrow(cohort) + value), ]$Pts <- count 
    }
      
      
    }
    #update the number of involved patients
    pts <- pts - discarded
    #check the number of patients
    if(pts >= (n_max - cohortsize + 1)){break}
    
    #after having assigned all the new pts we can take the decision for the next cohort --> based on all the data up to now (by the rrival of the new pts)
    results <- decision_EWOC(cohort, doses_info, time_arrival = last_arrival, target = target,  run = run, i_simulation = i_simulation)
    diagnostics <- rbind(diagnostics, results$diagnostics)
    info_probs <- results$info_probs
    next_dose <- results$choice
    
    if (is.na(next_dose)){
        n_early_stop <- n_early_stop + 1
        break
    }
    
    
    #check for SUfficient information (same as 3cohorts + 'STAY')
    if (next_dose == current_dose){
      Pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
      if(Pts$pts >= n_stop){
        n_stop_reached <- n_stop_reached + 1
        break
      }
    }
    
   
    #after the assignement of the new dose you shoudl check again for the bakcfill bcs you open those that are lower than the actual cohort dose (but do not open those with state 2 and state -1) --> if you add pts then you might have response so you have to open again. --> is it in the BOIN? IF i add pts at a dose then i do not check again for the new ppl if i can open again 
    
    if(current_dose <= doses_info$Dose[nrow(doses_info)]){
      index <- which(doses_info$Dose == current_dose)
      if (doses_info[index, ]$State == -1){
        parms_resp <- doses_info %>% filter(Dose == current_dose) %>% select(shape = Shape_resp, scale = Scale_resp)
        
          n_resp <- 0
          for (i in (pts-cohortsize):(pts - 1) ){
          resp <- rweibull(1, shape = parms_resp$shape, scale = parms_resp$scale)
          resp_i <- ifelse(resp < min((DLT_time + as.numeric(cohort[(i+1), ]$Time_arrival)), limit_time), 1, 0)
          n_resp <- n_resp + resp_i
        }
        
         #check cap for the backfill (if you can open but it has limit pts then put at 2)
          pts_check <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
          doses_info[index, ]$State <- ifelse(n_resp >= 1, ifelse(pts_check$pts < n_cap, 1, 2), -1)
      
      }   
    }
    
  #open the doses of backfill for safety (all those that are < next_dose and are not -1 or 2) and close all those for not-safety
    for (dose_state in 1:nrow(doses_info)){
      if (doses_info[dose_state, ]$Dose <= next_dose){
        if (doses_info[dose_state, ]$State == 0){
          doses_info[dose_state, ]$State <- 1
          #since we are ok with safety up to the MTD (we need to check for the number of pts? No since already checked within the backfill loop)
        }
      } else{
          if (doses_info[dose_state, ]$State == 1){
            doses_info[dose_state, ]$State <- 0 #for safety
          }
      }
    }
  

}

  #update the current dose based on the decision and the n_tot counter
  current_dose <- next_dose

}
  
  return(list('cohort' = cohort, 'diagnostics' = diagnostics, 'safety' = n_early_stop, 'sufficient' = n_stop_reached , 'doses_info' = doses_info))
}


#BLRM()
BFBLRM <- function(doses, true_pDLT, true_presp, cohortsize,  n_max, n_cap, n_stop, DLT_time, lambda, target = 0.3, i_simulation = 0){
  
  #compute the parms for each of the true_pDLT
  doses_info <- data.frame('Dose' = doses, 'Prob_DLT' = true_pDLT, 'Shape_DLT' = rep(NA, length(doses)), 'Scale_DLT' = rep(NA, length(doses)), 'Prob_resp' = true_presp, 'Shape_resp' = rep(NA, length(doses)), 'Scale_resp' = rep(NA, length(doses)))
  for(i in 1:nrow(doses_info)){
    
    parms_DLT <- weibull_parms(doses_info$Prob_DLT[i], DLT_time)
    doses_info$Shape_DLT[i] <- parms_DLT[1]
    doses_info$Scale_DLT[i] <- parms_DLT[2]
    
    parms_resp <- weibull_parms(doses_info$Prob_resp[i], DLT_time)
    doses_info$Shape_resp[i] <- parms_resp[1]
    doses_info$Scale_resp[i] <- parms_resp[2]
    
  }
  
  #add the status 
  doses_info$State <- rep(NA, length(doses)) #state of opening of the backfilling dose levels: NA at the beginning, 0 if closed, 1 if open, 2 if close since n_cap
  cohort <- data.frame('Run' = as.numeric(rep(NA, n_max)), 'Pts' = as.numeric(rep(NA, n_max)), 'Group' = rep(NA, n_max), 'Dose'= as.numeric(rep(NA, n_max)), 'Time_arrival' = as.numeric(rep(-1, n_max)), 'Time_DLT' = as.numeric(rep(-1, n_max)), 'DLT' = as.numeric(rep(-1, n_max)), 'Lag' = as.numeric(rep(-1, n_max)), 'Limit_time' = as.numeric(rep(-1, n_max)), 'Sum_times' = as.numeric(rep(-1, n_max)))
  results <- compute_BFBLRM_EWOC(cohort, doses_info, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda, target = target, i_simulation = i_simulation)
  return(results)
}

#---------------------------------------------------------------------------------------------------------------------

#Example
doses <- c(1.5, 2.5, 3.5, 4.5, 6, 7)
cohortsize <- 3
n_cap <- 12
n_stop <- 9
n_max <- 32
DLT_time <- 1 #1 month
lambda <- 6 


# Scenario 0
reference_dose <- 0.75
true_pDLT <- c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55, 0.6)
# Scenario 1
reference_dose <- 1.51
true_pDLT <- c(0.25, 0.41, 0.45, 0.49, 0.53, 0.57)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55, 0.6)
# Scenario 2
reference_dose <- 2.51
true_pDLT <- c(0.12, 0.25, 0.42, 0.49, 0.55, 0.62)
true_presp <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
# Scenario 3
reference_dose <- 3.51
true_pDLT <- c(0.04, 0.12, 0.25, 0.43, 0.63, 0.75)
true_presp <- c(0.1, 0.2, 0.3, 0.45, 0.58, 0.67)
# Scenario 4 
reference_dose <- 4.51
true_pDLT <- c(0.02, 0.06, 0.1, 0.25, 0.4, 0.5)
true_presp <- c(0.05, 0.1, 0.15, 0.3, 0.45, 0.55)
# Scenario 5
reference_dose <- 6.1
true_pDLT <- c(0.02, 0.05, 0.08, 0.11, 0.25, 0.39)
true_presp <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)
#scenario 6
reference_dose <- 7.1
true_pDLT <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.3)
true_presp <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)

#try new scenario a scenario 6 with higher probs of response at lower doses and see 
reference_dose <- 7.1
true_pDLT <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.3)
true_presp <- c(0.3, 0.35, 0.35, 0.4, 0.41, 0.44)

response_10 <- BFBLRM(doses, true_pDLT, true_presp, cohortsize,  n_max, n_cap, n_stop, DLT_time, lambda, i_simulation = 4)
