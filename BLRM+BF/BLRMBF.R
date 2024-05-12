library(dplyr)
library(tidyverse)
library(BOIN)
library(mvtnorm)
library(coda)
-----------------------------------------------------------------------------------------------------
## weibull_parms()

weibull_parms <- function(pDLT, DLT_time){
  #x1 = 0.5*(DLT_time) x2 = DLT_time
  shape <- (log(-log(1 - pDLT)) - log(-log(1 - 0.5*pDLT)))/(log(DLT_time) - log(0.5*DLT_time) )
  scale <- (0.5*DLT_time)/((-log(1 - 0.5*pDLT))^(1/shape))
  
  return(c(shape, scale))
}

## cohort_patient()

cohort_patient <- function(cohort, doses_info, current_dose, run, pts, previous_time, last_arrival_time = 0, last_lag = 0, i_simulation = 0){
  
  
  seed <- 1234 + run #+ i_simulation
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

## backfill_patients()

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

## decision_overtoxicity()
decision_overtoxicity <- function(Pts, DLT, target = 0.3, cutoff.eli = 0.95){
  
  boundaries <- get.boundary(target = target, ncohort = 10, cohortsize = cohortsize, cutoff.eli = cutoff.eli)$full_boundary_tab[4, ]
  
  if(DLT >= boundaries[Pts]){
    return(0)
  } else { return(1) }
}

## open_close()
open_close <- function(cohort, time_pts_backfill, dose, doses_info, target = 0.3){
  
  doses <- doses_info$Dose
        #open or closed?  
        
        #a: backfill and cohort at same dose --> of those whose Limit_time is less than the current time_arrival 
        DLT_tot_same_dose <- cohort %>% filter(Dose == dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(DLT = sum(DLT))
        pts_tot_same_dose <- cohort %>% filter(Dose == dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(Pts = n())
        
        #decide to keep open or to close 
        decision_overtoxicity_same_dose <- decision_overtoxicity(pts_tot_same_dose$Pts, DLT_tot_same_dose$DLT, target = target)
        
        #close the current dose of backfilling
        if(decision_overtoxicity_same_dose == 0) {
          #print('return 0')
          return(0)
          
        } else {
          return(1)
        }
}

## maximum_open()

maximum_open <- function(cohort, time_pts_backfill, current_dose, doses_info, target = 0.3){
  
      doses <- doses_info$Dose
  
      change_backfill <- F
      while (isFALSE(change_backfill)){
        
        #first condition: open state
        current_backfill_dose <- doses_info %>% filter(State == 1)  %>% filter(Dose < current_dose) %>% arrange(Dose) %>% tail(1)
        current_backfill_dose <- current_backfill_dose$Dose
  
        #valid dose?
        if (identical(current_backfill_dose, numeric(0))) {
          current_backfill_dose <- NA
          break
        }
        
        #open or closed?  Find the minimum dose open for backfill 
        response <- open_close(cohort, time_pts_backfill, current_backfill_dose, doses_info, target = target)
       
        #close the current dose of backfilling and all the upper doses 
        if(response == 0) {
          state <- doses_info$State
          state[which(doses == current_backfill_dose): length(state)] <- 0
          doses_info$State <- state
          
        } else {
          change_backfill <- T
        }
      }
      
      return(list('current_backfill_dose' = current_backfill_dose, 'doses_info' = doses_info))
}


## logposterior()
logposterior <- function(data, betas){
  
  #print the given values for beta0 and log(beta1)
  
  #define the prior parameters as in the BRLM+BF paper (6 doses)
  prior_mean <- c(log(0.5), 0)
  prior_var <- matrix(c(4, 0 ,0, 1), nrow = 2)
  
  #vector of probabilities based on the betas and on the doses with the model dependign of beta0 and beta1 (transform with exp)
  logit_p <- betas[1] + exp(betas[2])*log(data$Doses/reference_dose)
  p <- 1/(1 + exp(-logit_p)) 

  likelihood <- sum(dbinom(data$DLT, data$Pts, p, log = T))
  #prior for the values of beta0 and log(beta1)
  prior <- dmvnorm(betas, mean = prior_mean, sigma = prior_var, log = T)

  #multiply for the jacobian (exp(beta[2])) to obtain the posterior in terms of beta0 and log(beta1). SInce log then it is beta2
  return(sum(sum(likelihood, prior), betas[2]))
  
}

## MCMC()

MCMC <- function(cohort, doses_info, time_arrival, i_simulation, run, iterations = 10000, burnin = 1000){
  
  #set the seed according to the run and to the simulation number. Ok global variables  
  #set.seed((1234 + run + i_simulation))
  set.seed(1234 + run)
  
  #define the sd factor for the acceptance rate
  sd <- 2.3
  
  
  #take the values in a summary format
  dlt <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(freq = sum(DLT))
  pts <- cohort %>% filter(Limit_time <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(pts = n())
  doses <- cohort %>% filter(!is.na(Dose)) %>% reframe(dose = unique(Dose))
  
  data <- data.frame('DLT' = dlt$freq, 'Pts' = pts$pts, 'Doses' = doses$dose )
  
  #print(c('data', data))
  mean_current <- c(0, 1)
  beta0_start <- 1
  beta1_start <- 0
  sigma_start <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

  
  #print(c('parms', beta0_start, beta1_start, sigma_start))
  
  #implement MCMC calling the logpost function 
  
  #parameters 
  n_accept <- 0
  beta_parms <- matrix(nrow = iterations, ncol = 2)
  colnames(beta_parms) <- c('beta0', 'beta1')
  beta_parms[1, ] <- c(beta0_start, beta1_start)

  lambda_proposed <- log((2.38^2)/2)
  
  for( i in 2:iterations){
  
    beta_proposal <- rmvnorm(n = 1, mean = mean_current, sigma = exp(lambda_proposed)*sigma_current) #sampling of beta0 and log(beta1)
     
    log_now <- logposterior(data, beta_proposal)
    log_previous <- logposterior(data, beta_parms[(i - 1), ])
    log_ratio <- log_now - log_previous
    #print(c('logratio', log_ratio))

    #this shoudl not be entered
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
    
    #updating 
    
    #weights
    gamma <- 1/(1 + (i - 1))^delta 
    #update other parms 
    #lambda tuning on the logscale
    lambda_proposed <- gamma*(min(exp(log_ratio), 1) - 0.24) + lambda_proposed
    #sigma varcov --> add a small value not to get stuck?
    sigma_current <- sigma_current + gamma*((beta_parms[i, ] - mean_current)%*%t(beta_parms[i, ] - mean_current) - sigma_current ) + 0.001
    #mean vector
    mean_current <- mean_current + gamma*(beta_parms[i, ] - mean_current)
    
  }
  
  #trim the dataset for the burnin 
  beta_parms <- beta_parms[-c(1:burnin), ]
  
  #run and collect the diagnostics: Geweke and acceptance rate 
  diagnostics <- c(round(n_accept/(iterations - 1), digits = 2) * 100, geweke.diag(as.mcmc(beta_parms[, 1]))$z,  geweke.diag(as.mcmc(beta_parms[, 2]))$z )

  #transform the beta1 direclty 
  beta_parms[, 2] <- exp(beta_parms[, 2])
  
  return(list('betas' = beta_parms, 'diagnostics' = diagnostics))
  
}
 
## decision()
decision <- function(cohort, doses_info, time_arrival = 1000, run, target = 0.3, i_simulation = 0){
  
  #directly call the MCMC
  results <- MCMC(cohort, doses_info, time_arrival, i_simulation = i_simulation, run = run)
  prob <- data.frame(matrix(ncol = length(doses_info$Dose), nrow = 9000))
  
  #compute the probabilities (exp the beta_1 as you have samples of log(beta_1))
  for (i in 1:length(doses_info$Dose) ){
    logit_prob <- results$betas[, 1] + results$betas[, 2]*log(doses_info$Dose[i]/reference_dose)
    prob[, i] <- exp(logit_prob)/(1 + exp(logit_prob))
    #hist(prob[, i])
  }
  #compute the abs of the differences --> what if you have two? take the lowest possibly (that one coming from the lower side: conservative --> ask Alexandre or look for it on the net)  
  distances <- abs(prob - target)
  print(c('distances', distances))
  choice <- doses_info$Dose[which(distances == min(distances))]
  #return the lowest in case you have double choice 
  return(list('choice' = choice[1], 'diagnostics' = results$diagnostics))
  
}

## k_fold_skipping()
k_fold_skipping <- function(current_dose, next_dose, doses, k = 2){
  
  maximum_dose <- current_dose*3
  
  next_dose <- ifelse(next_dose <= maximum_dose, next_dose, maximum_dose)
  
  return(next_dose)
}

## hard_safety()

hard_safety <- function(cohort, time_arrival = 1000, target = 0.3){
  
  Pts <- cohort %>% filter(Time_DLT <= time_arrival) %>% filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(pts = n())
  DLT <- cohort %>% filter(Time_DLT <= time_arrival) %>%  filter(!is.na(Dose)) %>% group_by(Dose) %>% summarise(DLT = sum(DLT)) 
  
  result <- NA
  
  for ( dose in Pts$Dose ){
    
    pts_dose <- Pts %>% filter(Dose == dose) %>% select(pts)
    DLT_dose <- DLT %>% filter(Dose == dose) %>% select(DLT)
    
    if (decision_overtoxicity(pts_dose$pts, DLT_dose$DLT, target) == 0){
      result <- dose
      break
    }

  }
  
  return(result)
  
}

## check_n_cap()
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

## transform_data()
transform_data <- function(cohort){
  
  suppressWarnings(cohort2 <- transform(cohort, Run = as.numeric(Run), Pts = as.numeric(Pts) , Dose = as.numeric(Dose), Time_arrival = as.numeric(Time_arrival), Time_DLT = as.numeric(Time_DLT), DLT = as.numeric(DLT), Lag = as.numeric(Lag), Limit_time = as.numeric(Limit_time), Sum_times = as.numeric(Sum_times) ))
  
  return(cohort2)
  
}

## transform_data_doses()
transform_data_doses <- function(doses_info){
  
  suppressWarnings(cohort2 <- transform(doses_info, Dose = as.numeric(Dose), Shape_DLT = as.numeric(Shape_DLT) , Scale_DLT = as.numeric(Scale_DLT), Prob_resp = as.numeric(Prob_resp), Shape_resp = as.numeric(Shape_resp), Scale_resp = as.numeric(Scale_resp), State = as.numeric(State)))
  
  return(cohort2)
  
}

## compute_BRBLRM()

compute_BFBLRM <- function(cohort, doses_info, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda, target = 0.3 ,i_simulation = 0){
  
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

  print(c('use data', doses_info))
  
  while (pts <= n_max){
  
    run <- run + 1
   print(c('currentdose', current_dose))
    
    #no backfilling
    if (sum(is.na(doses_info$State)) != 0 ){
      print('No backfill open yet')
     
      print(c('current pts prior to cohort call', pts))
    update <- cohort_patient(cohort, doses_info, current_dose, run, pts, previous_time = limit_time)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    previous_time <- update$previous_time
    limit_time <- update$limit_time
    #check the number of pts
    print(pts)
    if(pts >= (n_max - cohortsize + 1)){break}
   
    #compute the decision and save the diagnostics 
    print('call decision')
    results <- decision(cohort, doses_info, target = target, run = run, i_simulation = i_simulation)
    diagnostics <- rbind(diagnostics, results$diagnostics)
    print(c('diagn',diagnostics))
    print(c('choice', results$choice))
    #deciison + safety rules at section 3.3
    next_dose <- results$choice
    #check if the next_dose is ok with the K-fold skipping Dose rule (this depends on the values of the doses. If you use indeces from 1 to 6 then you might have problems with this rule )
    next_dose <- k_fold_skipping(current_dose, next_dose, doses)
    
    #check the hard safety for all the doses and eliminate all those that are over the limit_dose (including this)
    limit_dose <- hard_safety(cohort, target = target) #this dose must be removed
    print(c('limit_dose', limit_dose))
    if (!is.na(limit_dose)){
      print('remove data')
      doses_info <- doses_info[-c(which(doses == limit_dose):nrow(doses_info)), ]
      doses <- doses_info$Dose 
     
        #check if the limit dose is not the minimum dose (apply the same limit of the BOIN alias 95%. This is different from the paper BLRM)
      if (limit_dose == min(doses)){
        print('Stop the trial')
        n_early_stop <- n_early_stop + 1
        break
      }
     
    }
    
    #check if the next_dose is still there available
    next_dose <- ifelse(next_dose <= doses_info$Dose[nrow(doses_info)], next_dose, doses_info$Dose[nrow(doses_info)])
    print(c('new next dose', next_dose))
    
    #check for SUfficient information (same as 3cohorts + 'STAY')
    if (next_dose == current_dose){
      Pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
      if(Pts$pts == n_stop){
        print('Maximum find')
        n_stop_reached <- n_stop_reached + 1
        break
      }
    }
    
    #if the next dose for the cohort is not the minimum dose then the backfill can start; if the current dose is the maximum one then the backfill option is not explored for the current dose level 
    backfill <- ifelse(next_dose != min(doses), ifelse(current_dose != max(doses), T, F), F)
    
          
    #possibility of backfilling
    if(backfill == T){
      print('possibility for backfill')
    
      parms_resp <- doses_info %>% filter(Dose == current_dose) %>% select(shape = Shape_resp, scale = Scale_resp)
    
      n_resp <- 0
    for (i in (pts-cohortsize):(pts - 1) ){
      resp <- rweibull(1, shape = parms_resp$shape, scale = parms_resp$scale)
      resp_i <- ifelse(resp < min((DLT_time + as.numeric(cohort[(i+1), ]$Time_arrival)), limit_time), 1, 0)
      n_resp <- n_resp + resp_i
    }
    
    if (n_resp >= 1){
      state <- ifelse(doses_info$Dose == current_dose, 1, ifelse(doses_info$Dose < current_dose, -1, 0))
      doses_info$State <- state
    }
}
  } else {
    
    #suppose the backfill is open so then we can work with this 
    
    print('backfill')#}
  
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
    
    
    #how to check if the backfill dose is still ok? we have the n_cap for the state 2 and the response for the state 1 but for closing the dose you have the state 0: what if we close a dose for overtoxicity but then we suspend only the backfill adn not the main cohort? then you'd use the values of the pts in the backfill for the final decision but you do not expose more people thna needed at the eventually bad situation --> but this would be based on the hard safety rule + i would do this for each of the pts in the backfill but for the pts in the cohort i'd stay with blocks of 3 --> either you close a dose for pts or for overtoxicity that you see at that moment but not for the de-escalation (that is based on the regression and here you do not have boundaries)
    discarded <- 0
    for (back_pts in temp_pts){
      print('backfill on')
      #actual time of arrival to be considered to check if the current dose for backfill is stll okay 
      time_pts_backfill <- cohort %>% filter(Pts == back_pts) %>% select(time = Time_arrival)
      update_backfill <- maximum_open(cohort, time_pts_backfill$time, current_dose, doses_info, target = target)
      current_maximum_backfill_dose <- update_backfill$current_backfill_dose
      doses_info <- transform_data_doses(update_backfill$doses_info)
      print(c('new patient at dose', current_maximum_backfill_dose))
      print(c('line 867', doses_info$State))
      
      #no backfill doses is open by arrival of the pts? then we discard the patient 
      if (is.na(current_maximum_backfill_dose)){
        print('No current_maximum')
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
    print(c('line 906', pts))
    if(pts >= (n_max - cohortsize + 1)){break}
    
    #after having assigned all the new pts we can take the decision for the next cohort --> based on all the data up to now (by the rrival of the new pts)
     print('call decision')
    results <- decision(cohort, doses_info, time_arrival = last_arrival, target = target,  run = run, i_simulation = i_simulation)
    diagnostics <- rbind(diagnostics, results$diagnostics)
    print(c('backfill diagn', diagnostics))
    next_dose <- results$choice
    print(c('backfill next dose', next_dose, 'current_dose', current_dose))
    #check if the next_dose is ok with the K-fold skipping Dose rule
    next_dose <- k_fold_skipping(current_dose, next_dose, doses)
    
    
    #check the hard safety for all the doses and eliminate all those that are over the limit_dose (including this)
    limit_dose <- hard_safety(cohort, target = target) #this dose must be removed
    print(c('backfill limit dose', limit_dose))
    if (!is.na(limit_dose)){
      print('removal')
      doses_info <- doses_info[-c(which(doses == limit_dose):nrow(doses_info)), ]
      doses <- doses_info$Dose 
     
        #check if the limit dose is not the minimum dose (apply the same limit of the BOIN alias 95%. This is different from the paper BLRM)
      if (limit_dose == min(doses)){
        print('Stop the trial')
        n_early_stop <- n_early_stop + 1
        break
      }
     
    }
    
    #check if the next_dose is still there available
    next_dose <- ifelse(next_dose <= doses_info$Dose[nrow(doses_info)], next_dose, doses_info$Dose[nrow(doses_info)])
    print(c('next dose check', next_dose))
    
    #check for SUfficient information (same as 3cohorts + 'STAY')
    if (next_dose == current_dose){
      Pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
      if(Pts$pts == n_stop){
        print('Maximum find')
        n_stop_reached <- n_stop_reached + 1
        break
      }
    }
    
    #after the assignement of the new dose you shoudl check again for the bakcfill bcs you open those that are lower than the actual cohort dose (but do not open those with state 2 and state -1) --> if you add pts then you might have response so you have to open again. --> is it in the BOIN? IF i add pts at a dose then i do not check again for the new ppl if i can open again 
    
    
    
    #check for opening of previous not-responsive in case the current dose has not been eliminated 
    if(current_dose <= doses_info$Dose[nrow(doses_info)]){
      print('ok try opening')
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
      print(c('current dose state', doses_info[dose_state, ]$State, 'at dose', dose_state))
      if (doses_info[dose_state, ]$Dose <= next_dose){
        if (doses_info[dose_state, ]$State == 0){
          print('open')
          doses_info[dose_state, ]$State <- 1
          #since we are ok with safety up to the MTD (we need to check for the number of pts? No since already checked within the backfill loop)
        }
      } else{
          if (doses_info[dose_state, ]$State == 1){
            print('close')
            doses_info[dose_state, ]$State <- 0 #for safety
          }
      }
    }

}

  #update the current dose based on the decision and the n_tot counter
  current_dose <- next_dose

}
  
  return(list('cohort' = cohort, 'diagnostics' = diagnostics))
}

## BFBLRM()

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

  results <- compute_BFBLRM(cohort, doses_info, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda, target = target, i_simulation = i_simulation)


  return(results)
}

---------------------------------------------------------------------------------------------------------------------------------------
# Example

doses <- seq(1, 6, 1)
cohortsize <- 3
n_cap <- 12
n_stop <- 9
n_max <- 12
DLT_time <- 1 #1 month
lambda <- 6 
#scenario 2
true_pDLT <- c(0.12, 0.25, 0.42, 0.49, 0.55, 0.62)
true_presp <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

BFBLRM(doses, true_pDLT, true_presp, cohortsize,  n_max, n_cap, n_stop, DLT_time, lambda)
