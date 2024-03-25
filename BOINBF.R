# load libraries 
library(dplyr)
library(tidyverse)
library(BOIN)

-------------------------------------------------------------------------------------------------------------------------

# Baseline functions 

## weibull_parms()
weibull_parms <- function(pDLT, DLT_time){
  #x1 = 0.5*(DLT_time) x2 = DLT_time
  shape <- (log(-log(1 - pDLT)) - log(-log(1 - 0.5*pDLT)))/(log(DLT_time) - log(0.5*DLT_time) )
  scale <- (0.5*DLT_time)/((-log(1 - 0.5*pDLT))^(1/shape))
  
  return(c(shape, scale))
}


## cohort_patient()
cohort_patient <- function(cohort, doses_info, current_dose, run, 
                           pts, previous_time, last_arrival_time = 0, last_lag = 0){
  
  # seed <- 1234 + run
  # set.seed(seed)
  
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

## decision()
decision <- function(new_reference, current_cohortsize, DLT){
  
  bounds <- new_reference %>% filter(Pts == current_cohortsize)
  change <- NULL
  
if (DLT >= bounds$Eliminate){
  change = 4
} else if (DLT <= bounds$Escalate) { change = 1} else if (DLT >= bounds$Descalate) { change = 3 } else {change = 2}
  
  return(change)
  
}

## backfill_patients()
backfill_patients <- function(cohort, run, pts, previous_time, limit_time, n_max){
    
    # seed <- 555+run
    # set.seed(seed)
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

## open_close()

open_close <- function(cohort, time_pts_backfill, dose, new_reference, doses_info){
  
  doses <- doses_info$Dose
        #open or closed?  
        
        #a: backfill and cohort at same dose --> of those whose Limit_time is less than the current time_arrival 
        DLT_tot_same_dose <- cohort %>% filter(Dose == dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(DLT = sum(DLT))
        pts_tot_same_dose <- cohort %>% filter(Dose == dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(Pts = n())
        
        #b: backfill and cohort at (same dose + 1) 
        upper_dose <- doses[(which(doses == dose)) + 1]
        DLT_tot_upper_dose <- cohort %>% filter(Dose == upper_dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(DLT = sum(DLT))
        pts_tot_upper_dose <- cohort %>% filter(Dose == upper_dose) %>% filter(Limit_time <= time_pts_backfill) %>% summarise(Pts = n())
        
        #do open the dose if there is no pts available for decision : it is not altering the level of the backfill if no info available
        if (pts_tot_same_dose == 0 || pts_tot_upper_dose == 0) {
          
          return(1)
        }
        
        #decide to keep open or to close 
        decision_same_dose <- decision(new_reference, pts_tot_same_dose$Pts, DLT_tot_same_dose$DLT)
        decision_pooled_dose <- decision(new_reference, ( pts_tot_upper_dose$Pts + pts_tot_same_dose$Pts), (DLT_tot_upper_dose$DLT + DLT_tot_same_dose$DLT) )
       
        #close the current dose of backfilling
        if((decision_same_dose == 3 || decision_same_dose == 4) && (decision_pooled_dose == 3 || decision_pooled_dose == 4) ) {
          return(0)
          
        } else {
          return(1)
        }
}

## maximum_open()

maximum_open <- function(cohort, time_pts_backfill, current_dose, doses_info, new_reference){
  
      doses <- doses_info$Dose
  
      change_backfill <- F
      while (isFALSE(change_backfill)){
        
        #first condition: open state
        current_backfill_dose <- doses_info %>% filter(State == 1)  %>% filter(Dose < current_dose) %>% arrange(Dose) %>% tail(1)
        current_backfill_dose <- current_backfill_dose$Dose
      
        #valid dose?
        if (identical(current_backfill_dose, numeric(0))) {
          #1-print('No length')
          current_backfill_dose <- NA
          break
        }
        
        #open or closed?  Find the minimum dose open for backfill 
        response <- open_close(cohort, time_pts_backfill, current_backfill_dose, new_reference, doses_info)

        #close the current dose of backfilling and all the upper doses 
        if(response == 0) {
          #1-print('dose closed')
          state <- doses_info$State
          state[which(doses == current_backfill_dose): length(state)] <- 0
          doses_info$State <- state
          
        } else {
          change_backfill <- T
        }
      }
      
      return(list('current_backfill_dose' = current_backfill_dose, 'doses_info' = doses_info))
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


## merged_decision()

merged_decision <- function(new_reference, cohort, dose_cohort, dose_backfill, last_arrival, doses_info){
  
  conflict_table <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 ), nrow = 4, ncol = 4, byrow = T)
  
  doses <- doses_info$Dose
  
  if(!identical(dose_backfill, numeric(0))){
    backfill_pts <- cohort %>% filter(Dose == dose_backfill) %>% filter(Limit_time <= last_arrival) %>% summarise(pts = n())
    backfill_DLT <- cohort %>% filter(Dose == dose_backfill) %>% filter(Limit_time <= last_arrival) %>% summarise(DLT = sum(DLT) )
    
    change_backfill <- decision(new_reference, backfill_pts$pts, backfill_DLT$DLT)
  } else { 
    change_backfill <- NA
    }

  cohort_pts <- cohort %>% filter(Dose == dose_cohort) %>% summarise(pts = n())
  cohort_DLT <- cohort %>% filter(Dose == dose_cohort) %>% summarise(DLT = sum(DLT) ) 
  
  change_cohort <- decision(new_reference, cohort_pts$pts, cohort_DLT$DLT)
  
  if(is.na(change_backfill) || conflict_table[change_backfill, change_cohort] != 0){
    if (is.na(change_backfill) ){print('no backfill')} else if(conflict_table[change_backfill, change_cohort] != 0){print('no conflict')} 
     if (change_cohort == 4){
      safe_dose <- doses[(which(doses == dose_cohort)) - 1]
      doses <- doses[-c(which(doses == dose_cohort):length(doses))]
      if (length(doses) == 0){
        safe_dose <- NA
      }
    } else if (change_cohort == 1){
      
      if(dose_cohort != doses[length(doses)]){
        safe_dose <- doses[(which(doses == dose_cohort) + 1)]
      } else { safe_dose <- dose_cohort }
    } else if (change_cohort == 3){
      if (dose_cohort != doses[1]){
        safe_dose <- doses[(which(doses == dose_cohort) - 1)]
      } else {safe_dose <- dose_cohort}
    } else {
      safe_dose <- dose_cohort
    }
    
    return(list('safe_dose' = safe_dose, 'change' = change_cohort))
  } else {
    
    #take all the doses from the actual backfilled one to the cohort
    all_doses <- doses[which(doses == dose_backfill):(which(doses == dose_cohort) - 1)]
    all_doses <- sort(all_doses, decreasing = T) #analise the doses from the highest
    changes <- c()
    
    safe_dose <- NA
    
    for (try_dose in all_doses){
      #1-print(c('try doses', try_dose))
      if(try_dose != dose_backfill){
      
        cumulative_pts <- cohort %>% filter(Dose >= dose_backfill) %>% filter(Dose <= try_dose) %>% summarise(pts = n())
        cumulative_DLT <- cohort %>% filter(Dose >= dose_backfill) %>% filter(Dose <= try_dose) %>% summarise(DLT = sum(DLT))
      
        cumulative_change <- decision(new_reference, cumulative_pts$pts, cumulative_DLT$DLT)
      } else { cumulative_change <- change_backfill }
      
      if (cumulative_change == 2){
        safe_dose <- try_dose
        break 
      } else if (cumulative_change == 1){
        if (try_dose != max(doses)){
          safe_dose <- doses[which(doses == try_dose) + 1]
        } else { 
          safe_dose <- try_dose
          break
          }
    
      } else { next }#if 4 or 3 we pass and let the try_dose change and go through the if tests
      
    } 
    
    if (is.na(safe_dose)){
      if (dose_backfill != min(doses)){
        safe_dose <- doses[(which(doses == dose_backfill) - 1)]
      } else {
        safe_dose <- dose_backfill
      }
      
    }
    
    
    return(list('safe_dose' = safe_dose, 'change' = cumulative_change))
    
  }
  
}

## transform_data()

transform_data <- function(cohort){
  
  suppressWarnings(cohort2 <- transform(cohort, Run = as.numeric(Run), Pts = as.numeric(Pts) , Dose = as.numeric(Dose), Time_arrival = as.numeric(Time_arrival), Time_DLT = as.numeric(Time_DLT), DLT = as.numeric(DLT), Lag = as.numeric(Lag), Limit_time = as.numeric(Limit_time), Sum_times = as.numeric(Sum_times) ))
  
  return(cohort2)
  
}

## compute_BFBOIN()

compute_BFBOIN <- function(cohort, doses_info, new_reference, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda){
  
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
    
  #initiate the while loop that should be stopped asap we reach: 1) n_max, 2) overtoxicity at minimal dose, 3) n_stop + stay decision 

  
  while (pts <= n_max){
    
    run <- run + 1
    
    #no backfilling
    if (sum(is.na(doses_info$State)) != 0 ){
    
    backfill <- T 
    update <- cohort_patient(cohort, doses_info, current_dose, run, pts, previous_time = limit_time)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    previous_time <- update$previous_time
    limit_time <- update$limit_time
    if(pts > n_max || pts > (n_max - cohortsize)){break}
    
    
    #compute the decision 
    #a) total DLT at current dose
    total_cohort_DLT <- cohort %>% filter(Dose == current_dose) %>% filter(Group == 'C') %>% summarise(DLT = sum(DLT))
    
    #b) pts at current dose 
    total_cohort_pts <- cohort %>% filter(Dose == current_dose) %>% filter(Group == 'C') %>% summarise(pts = n())
    
    #c) compute the decision and save it (check purposes)
    change <- decision(new_reference, total_cohort_pts$pts, total_cohort_DLT$DLT)
    print(c('change', change))
    
    #d) evaluate the decision 
    if (change == 4){
      index <- (which(doses == current_dose))
      next_dose <- doses[(which(doses == current_dose)) - 1]
      doses_info <- doses_info[-c(index:nrow(doses_info)), ] #eliminate all the upper doses
      doses <- doses_info$Dose
      state <- ifelse(doses_info$Dose == current_dose, 0, doses_info$State)
      doses_info$State <- state
      backfill <- F
      #doses <- doses[-c(which(doses == current_dose):length(doses))]
      if (length(doses) == 0){
        n_early_stop <- n_early_stop + 1
        break
      }
    } else if (change == 1){
      
      if(current_dose != doses[length(doses)]){
        next_dose <- doses[(which(doses == current_dose) + 1)]
      } else { 
        next_dose <- current_dose
        #check the stop 
          if(total_cohort_pts == n_stop || total_cohort_pts > n_stop){
          n_stop_reached <- n_stop_reached + 1
          break 
         }
       }
    } else if (change == 3){
      state <- ifelse(doses_info$Dose == current_dose, 0, doses_info$State)
        doses_info$State <- state
        backfill <- F
      if (current_dose != doses[1]){
        next_dose <- doses[(which(doses == current_dose) - 1)]
      } else {
        #check the stop 
        next_dose <- current_dose
          if(total_cohort_pts == n_stop || total_cohort_pts > n_stop){
          n_stop_reached <- n_stop_reached + 1
          break 
        }
      }
    } else {
      next_dose <- current_dose
      if (current_dose == min(doses)){
        backfill <- F
      } 
      #check the stop
      if(total_cohort_pts == n_stop || total_cohort_pts > n_stop){
        n_stop_reached <- n_stop_reached + 1
        break 
      }
    }

    #maximum dose?
    if (current_dose == max(doses)){
      backfill <- F
    }
          
    #possibility of backfilling
    if(isTRUE(backfill)){
    
    #strategy 1: binomial 
      
    # p <- doses_info %>% filter(Dose == current_dose) %>% select(Prob_resp)
    # n_resp <- rbinom(1, n = cohortsize, p = p$Prob_resp)
    # 
    # doses_info$State[(which(doses == current_dose))] <- ifelse(n_resp >= 1, 1, 0)
    #  
    
    #strategy 2: use a weibull with DLT_time as maximum waiting 
    #a): simulate the response for each of the patient --> if it happens during his/her DLT_time 
    #b): if the response is after min(DLT_time + time_arrival, limit_time) then set it equal to 0 as you'd stop following the patient 
    
      parms_resp <- doses_info %>% filter(Dose == current_dose) %>% select(shape = Shape_resp, scale = Scale_resp)
    
      n_resp <- 0
    for (i in (pts-cohortsize):(pts - 1) ){
      resp <- rweibull(1, shape = parms_resp$shape, scale = parms_resp$scale)
      resp_i <- ifelse(resp < min((DLT_time + as.numeric(cohort[(i+1), ]$Time_arrival)), limit_time), 1, 0)
      n_resp <- n_resp + resp_i
    }
    
    index <- which(doses == current_dose)
    if (n_resp >= 1){
      state <- ifelse(doses_info$Dose == current_dose, 1, ifelse(doses_info$Dose < current_dose, -1, 0))
      doses_info$State <- state
    }
}
  } else {
    
  
    #vector of backfill has been explored already (no NA)
    
    # cohort pts --> cohortsize added
    update <- cohort_patient(cohort, doses_info, current_dose, run, pts, previous_time = limit_time, last_arrival_time = last_arrival, last_lag = last_lag)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    #previous_time for backfill adn limit_time for cohort
    previous_time <- update$previous_time
    limit_time <- update$limit_time
    #check the number of pts 
    #print(pts)
    if(pts > n_max){break}


    #compute the limit_time across the pts of the current run 
    # limit_time <- cohort %>% filter(Group == 'C') %>% filter(Run == run) %>% summarise(limit_time = max(Limit_time))
    # limit_time <- as.numeric(limit_time$limit_time)
    
    #simulate arrival times 
    update <- backfill_patients(cohort, run, pts, previous_time = previous_time, limit_time = limit_time, n_max)
    cohort <- transform_data(update$cohort)
    pts <- update$pts
    #previous_time <- update$previous_time
    #save last arrivals
    last_arrival <- update$last_arrival
    last_lag <- update$last_lag
    
    
    #decide the dose level for each of the arrived patients 
    temp_pts <- cohort %>% filter(Run == run) %>% filter(Group == 'B') %>% select(pts = Pts)
    temp_pts <- temp_pts$pts
    
    discarded <- 0
    for (back_pts in temp_pts){
      #actual time of arrival to be considered
      time_pts_backfill <- cohort %>% filter(Pts == back_pts) %>% select(time = Time_arrival)
      update_backfill <- maximum_open(cohort, time_pts_backfill$time, current_dose, doses_info, new_reference)
      current_maximum_backfill_dose <- update_backfill$current_backfill_dose
      doses_info <- update_backfill$doses_info
      
      #no backfill doses is open by arrival of the pts? Use min(doses) as reference to start the filtering from 
      if (is.na(current_maximum_backfill_dose)){
        higher_doses <- doses_info %>% filter(State != -1) %>% filter(State != 2) %>% filter(Dose >= min(doses)) %>% filter(Dose < current_dose) %>% select(Dose)
      higher_doses <- higher_doses$Dose
      } else {
        higher_doses <- doses_info %>% filter(State != -1) %>% filter(State != 2) %>% filter(Dose > current_maximum_backfill_dose) %>% filter(Dose < current_dose) %>% select(Dose)
      higher_doses <- higher_doses$Dose
      }
      
      #select doses that can be opened again 
      if (!(identical(higher_doses, numeric(0))) && sum(is.na(higher_doses)) == 0) {
        change <- T
        i <- 1
        while(i <= length(higher_doses) ){
          possible_dose_state <- open_close(cohort, time_pts_backfill$time, higher_doses[i], new_reference, doses_info)
          if (possible_dose_state == 0){
            #keep the dose closed
            break #it means that we can not open this and the upper ones 
          } else {
            # re-open the dose. the round's highest one is taken as the backfilling dose
            state <- doses_info$State
            state[which(doses == higher_doses[i])] <- 1
            doses_info$State <- state
            current_backfill_dose <- higher_doses[i]  
          }
           i <- i +1 
        }
      } else {
        current_backfill_dose <- current_maximum_backfill_dose
      }
      
      if (is.na(current_backfill_dose)){
        discarded <- discarded + 1 
        next
      }
      
      #assign the dose and compute the time_DLT
      parms_backfill_DLT <- doses_info %>% filter(Dose == current_backfill_dose) %>% select(shape = Shape_DLT, scale = Scale_DLT)
      time_DLT <- round(rweibull(1, shape = parms_backfill_DLT$shape, scale = parms_backfill_DLT$scale), 4)
      cohort[back_pts, ]$Dose <- current_backfill_dose
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
    if(pts > n_max || pts > (n_max - cohortsize)){break}
    
    #after having assigned all the new pts we can take the decision for the next cohort --> based on the highest backfill used 
    
    dose_backfill <- cohort %>% filter(Run == run) %>% filter(Group == 'B') %>% arrange(Dose) %>% tail(1)
    next_dose <- merged_decision(new_reference, cohort, current_dose, dose_backfill$Dose, last_arrival = last_arrival, doses_info) 
    
    #evaluation of the results for the cohort 
     
    if(is.na(next_dose$safe_dose)){
      n_early_stop <- n_early_stop + 1
      break
    } else if (next_dose$change == 2){
        #stay at the same dose: n_stop reached?
        total_cohort_pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
        if(total_cohort_pts$pts >= n_stop){
          n_stop_reached <- n_stop_reached + 1
          break 
        }
    } else if (change == 3 && current_dose == min(doses_info$Dose)){
       
        total_cohort_pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
        if(total_cohort_pts$pts >= n_stop){
          n_stop_reached <- n_stop_reached + 1
          break 
        }
    } else if (change == 1 && current_dose == max(doses_info$Dose)){
      
        total_cohort_pts <- cohort %>% filter(Dose == current_dose) %>% summarise(pts = n())
        if(total_cohort_pts$pts >= n_stop){
          n_stop_reached <- n_stop_reached + 1
          break 
        }
      
    } else if (change == 4){
      doses_info <- subset(doses_info, Dose <= next_dose$safe_dose)
      doses <- doses_info$Dose
    }    
     
  #check for the basckfill   
     
  if ((next_dose$change == 3 || next_dose$change == 4) && sum(doses_info$State) != 0 ) {
      #close upper 
      state <- ifelse(doses_info$Dose >= current_dose, 0, doses_info$State)
      doses_info$State <- state
    } 

        
        
    #keep just the dose 
    next_dose <- next_dose$safe_dose
    
}

  #update the current dose based on the decision and the n_tot counter
  current_dose <- next_dose

}
  
  return(list('cohort' = cohort, 'doses_infr' = doses_info))
}

## BFBOIN()

BFBOIN <- function(doses, true_pDLT, true_presp, cohortsize,
                   n_max, n_cap, n_stop, DLT_time, lambda, target = 0.3){
  
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
 
  #get the reference table
  reference_general <- get.boundary(target = 0.3, ncohort = 10, cohortsize = 3)$full_boundary_tab
  #manage to easier decision() handling
  new_reference <- as.data.frame(t(reference_general))
  colnames(new_reference) <- c('Pts', 'Escalate','Descalate','Eliminate')

  cohort <- data.frame('Run' = as.numeric(rep(NA, n_max)), 'Pts' = as.numeric(rep(NA, n_max)), 'Group' = rep(NA, n_max), 'Dose'= as.numeric(rep(NA, n_max)), 'Time_arrival' = as.numeric(rep(-1, n_max)), 'Time_DLT' = as.numeric(rep(-1, n_max)), 'DLT' = as.numeric(rep(-1, n_max)), 'Lag' = as.numeric(rep(-1, n_max)), 'Limit_time' = as.numeric(rep(-1, n_max)), 'Sum_times' = as.numeric(rep(-1, n_max)))

  results <- compute_BFBOIN(cohort, doses_info, new_reference, cohortsize, n_max, n_cap, n_stop, DLT_time, lambda)


  return(results)
}

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Example 

doses <- seq(1, 5, 1)
cohortsize <- 3
n_cap <- 12
n_stop <- 9
n_max <- 32
DLT_time <- 1
lambda <- 3

# Scenario 0
true_pDLT <- c(0.40, 0.45, 0.50, 0.55, 0.60)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55)
# Scenario 1
true_pDLT <- c(0.25, 0.41, 0.45, 0.49, 0.53)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55)
# Scenario 2
true_pDLT <- c(0.12, 0.25, 0.42, 0.49, 0.55)
true_presp <- c(0.2, 0.3, 0.4, 0.5, 0.6)
# Scenario 3
true_pDLT <- c(0.04, 0.12, 0.25, 0.43, 0.63)
true_presp <- c(0.1, 0.2, 0.3, 0.45, 0.58)
# Scenario 4 
true_pDLT <- c(0.02, 0.06, 0.1, 0.25, 0.4)
true_presp <- c(0.05, 0.1, 0.15, 0.3, 0.45)
# Scenario 5
true_pDLT <- c(0.02, 0.05, 0.08, 0.11, 0.25)
true_presp <- c(0.05, 0.1, 0.15, 0.2, 0.3)
# Scenario 6
true_pDLT <- c(0.12, 0.25, 0.42, 0.49, 0.55)
true_presp <- c(0.3, 0.35, 0.36, 0.36, 0.36)
# Scenario 7
true_pDLT <- c(0.04, 0.12, 0.25, 0.43, 0.63)
true_presp <- c(0.15, 0.3, 0.35, 0.36, 0.36)
# Scenario 8
true_pDLT <- c(0.02, 0.06, 0.1, 0.25, 0.4)
true_presp <- c(0.10, 0.2, 0.3, 0.35, 0.35)
# Scenario 9
true_pDLT <- c(0.02, 0.05, 0.08, 0.11, 0.25)
true_presp <- c(0.10, 0.15, 0.2, 0.3, 0.35)
# Scenario 10
true_pDLT <- c(0.04, 0.12, 0.25, 0.43, 0.63)
true_presp <- c(0.3, 0.32, 0.35, 0.36, 0.36)
# Scenario 11
true_pDLT <- c(0.02, 0.06, 0.1, 0.25, 0.4)
true_presp <- c(0.1, 0.3, 0.32, 0.35, 0.36)
# Scenario 12
true_pDLT <- c(0.02, 0.05, 0.08, 0.11, 0.25)
true_presp <- c(0.1, 0.15, 0.3, 0.32, 0.35)

results <- BFBOIN(doses, true_pDLT, true_presp, cohortsize,  n_max, n_cap, n_stop, DLT_time, lambda)

## MTD finding with select.MTD from BOIN

new_results <- transform_data(results$cohort) 
y <- new_results %>% group_by(Dose) %>% summarise(y = sum(Pts))
tox <- new_results %>% group_by(Dose) %>% summarise(tox = sum(DLT))

y <- subset(y, !is.na(Dose))
tox <- subset(tox, !is.na(Dose))

select.mtd(target = 0.3, npts = y$y, ntox = tox$tox)

