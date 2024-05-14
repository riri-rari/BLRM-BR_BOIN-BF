# Introduction 

These codes are aimed to implement the BF+BLRM algorithm. 
Both codes are following the model in Barnet et al., but the BLRMBF.R code uses a the mean of the posterior probability witht the hard safety and k-fold-skipping rules, while the BLRMBF_EWCO.R code uses the interval approach with EWOC in Neuenschwander (up to now with some modifications). 

# BF+BLRM 

## Baseline functions 

### weibull_parms()

weibull_parms() is a function that takes in input the probability of an event (as DLT at a specific dose, ie pDLT) and the time within which the event might happen (DLT_time, that is the time of patient follow-up). The function assumes that the DLT_Time is the pDLT quantile of a Weibull function and that $0.5*DLT_time$ is the $0.5*pDLT$ quantile of a Weibull distribution. From the quantiles, the function computes the values of the Shape and Scale parameters (ref: https://www.johndcook.com/quantiles_parameters.pdf).

### cohort_patient()

The cohort_patient function is devoted to the computation of the values for the cohort patients. It takes in input  the cohort dataset (cohort), the info about the doses (doses\_info), the current dose value (current\_dose), the current run (run), the current patient count (pts), the current previous time (previous_time), the last arrival time (last\_arrival\_time, default = 0 otherwise eqaul to the time of the person who stopped the backfilling. See backfill\_patients()) and the last lag (last\_lag, default = 0 otherwise set to the last lag of the person who stopped the backfilling. See backfill\_patients()). The function returns the cohort dataset updated with the patients info, the updated patients count and the updated previous time.

The function retrieves the parameters of the Weibull distribution at the current dose from the doses\_info dataset. Then, it cycles over a cohortsize of patients: if the patient is the first one of the trial, then his/her time\_arrival and lag\_arrival are not computed (set to 0) and the previous time is kept to the original value (it should be passed as 0); if the patient is not the first one, then his/her lag of arrival is simulated from the Exp(lambda), the time\_arrival is computed as the sum of the previous\_time (that represents the time of arrival of the previous patient. The input value of previous\_time is the limit\_time of the previous cohort) summed to the lag\_arrival and the previous\_time is then updated to the current patient time\_arrival. If a value for the last\_arrival\_time is given and the call of the function is 1 (so it is the first patient at the current cohort), the time\_arrival of the patient is set equal to the given arrival time and the lag is set equal to the given lag time (this is thought to reusing the patient that has broken the bakfilling process. See backfill\_patients()).

Regardless of the patient's number, the time for the DLT is computed (time\_DLT), the DLT is computed and the limit_time is computed (this represents the point in time when the patient is stopped to be followed-up).

The computed information are returned as a list. A summary of the time values is in Appendix 2

### backfill_patients()

The backfill\_patients() function is aimed at the simulation of all those patients that should enter a backfill arm as they arrive between the last cohort patient and the limit time. It returns an updated version of the cohort dataset with the backfillable patients (info of their dose and pDLT-related characteristics are set to NA) and the time that has make the while loop ceasing (the over-limit-time-patient's arrival time). 

Time arrival constrain is less-equal to the limit_time. 
Patients are generated up to the n\_max limit if all the patients arrives within the limit time.

THe function resturns the updated cohort, the updated number of patients, the updated previous time (last time of arrival of a backfilled patient), the last time arrival (last\_arrival) and the last lag (last\_lag) corresponsing to the arrival time and the lag time of the patient that interrupted the backfill.

### decision_overtoxicity()

The decision_overtoxicity() function is aimed to compute the $P(\pi_{posterior} > 0.3) \ge 0.95$ (with target DLTrate = 0.3). As the rule is the same as that implemented in the BOIN, the same limit for the number of DLT at that dose is applied. The function returns 0 if the dose is overly-toxic and it returns 1 if the dose is okay. 

The function takes in input the number of Pts at that dose (Pts), the number of DLT at that dose (DLT), the target toxicity (target, with default 0.3) and the cutoff for the probability (cutoff.eli, by default set to 0.95).
Notes for me: It also uses the global variable cohortsize. 

### check_individual_posterior()
The check_individual_posterior() is a function aimed to compute the posterior probbailty of a probability rv (conjugacy between Beta and Binomial). It computes then the probability of being over the toxicity limit (Neuen 0.35 >= 0.25). The function is used to check, by the arrival of new backfiling patient, the goodness of the posterior probability for the cohort and for the backfilling arms. 
The function returns a boolean variable that is T if the posterior is fine and it is F if the posterior is not fine. 

### open_cose()

The open_close() function is used to decide whether a specific dose should be open or closed for backfilling. This is done in in two steps: 
- a) collection of the total patients and DLT events at the same dose level. Only patients with (time_arrival + DLT_time) $\le$ current backfilled patient's arrival time are included (the only info that we have at this time point)
- b) computation of the decision for overtoxicity (decision_overtoxicity() function) and evaluation of the results. If the dose is closed for backfilling, all the doses greater than this are closed too (see maximum_open()).

If a dose should be closed then a 0 is returned, while if it should be open, function returns 1.

The function takes in input the cohort dataset, the time of the current patient of backfll (time\_pts\_backfill), the dose that is currently under investigation (dose), the reference dataset for the boundaries (new\_reference. See BFBOIN()) and the dataset of info about the doses (doses\_info).

The function returns ether 1 or 0, as specified above.

### open_close_EWOC()

The open_cloe_EWOC() function is aimed to find the optimal dose for backfilling. It exploits the MCMC to find the best dose matching with the target dose with the constraint of being lower than the current cohort dose (limit_dose set to the current dose index in the call of evaluate_posterior() function). Being an MCMC run then the function returns the diagnostics an the matrix of posterior probabilities information at each dose level. It returns also the chosen backfill dose.

### maximum_open()

The maximum_open() function serves to find the maximum dose open for the backfill. The function loops over the doses (each time it takes the highest dose open for backfill and lower than the current cohort dose) until a dose can be kept open (change\_backfill switches from F to T). The function returns the found backfill dose and the updated doses\_info dataset.

The function is deemed to check if the previously highest dose open for backfill is still okay and to eventually close it. Since the function loops until a change of dose is not found, all those doses that were open for backfill prior to the arrival of the current patient and that by the arrival of the patient should be closed, are closed. The maximum dose open is then found. 
This takes into consideration all those patients who has experienced a DLT between the previous patient and the current patient arrival time, changing the evaluation of the doses for backfilling. In the case for which no more open doses are left for backfilling (current\_dose\_backfilling is an empty object), the while loop is broken and the function returns a NA value. In the case that the current dose is available, the state of it is checked by means of the open\_close() function. 
If this function returns 1 then the dose is left open, while if it returns 0 the dose is closed as well as all the doses that are higher than the current one. 

In the case that a dose is closed, then the cycle is entered again. If the dose is, instead, the right one (open\_close returns 1), the cycle is exit (change = T) and the current dose is returned. 

The function takes in input the cohort dataset, the time of the current patient of backfll (time\_pts\_backfill), the dose that is currently under investigation (current\_dose), the dataset of info about the doses (doses\_info) and the reference dataset for the boundaries (new\_reference. See BFBOIN()) 

### maximum_open_EWOC()

The maximum_open_EWOC() function is aimed at finding the optimal dose for the backfilling arm by the arrival of the next patient. The function firslty finds the highest dose open for backfill and then checks if by the arrival of the current patients (time_pts_backfill) the dose level is still okay or not. To make the decision, the cohort arm and the current backfill dose arm are checked for their posterior probability of having P(x >= 0.35) >= 0.25 by means of the check_individual_posterior() function. If both the arms are okay, then the current backfill dose is used for the current patient and returned. If either of the arms is not okay then the open_close_EWOC() function is called and the optimal dose for the backfill is found by measn of MCMC run. The doses_info dataset is then updated such that all the doses over the current backfill one are closed for backfilling and all those avilable for being open (not -1 and not 2, see Appendix()) are open. In this scenario the chosen backfill dose, the diagnostics, the probabilities info and the updated version of the doses_info dataset are returned. 


### logposterior()

The logposterior() function is aimed to compute the logarithm of the posterior distribution for the $\beta$ parameters that are used inthe logistic model as specified by the following model: 

$Logit(\pi_i) = \beta_0 + \beta_1log(\frac{d_i}{d_{ref}}) $, with $d_i$ the dose of subject i. 

The prior that is used is the following:

$(\beta_0, \log(\beta_1) ) \sim N((\log(0.5), 0), (4, 0, 0, 1))$ (ref BLRM+BF paper. to be changed). 

The function takes in input the data (data dataframe with summary information about the DLT and the PTs for each of the used doses) and the vector of betas that are sampled inside the MCMC loop (betas, see MCMC()). 

Assumptions used in the formula: 
- independent obsrvations (ok different patients) within same experiemtn (same conditions) --> so i can write Lik_2 * Lik_1 * prior to have the Posterior_2
- exchangeability (no matter the order given the parameters beta) --> so i can use all info at once

### MCMC()

The MCMC() function is aimed to compute the MCMC posterior densities estimates for the $\beta_0$ and $\log(\beta_1)$ parameters. It does so applying the Metropolis-Hastings algorithm. 

The function computed the DLT and the Pts for all the doses that have been explored so far (only patients with Limit_time below the arrival time are included) and it stores them in the data dataframe. Within the MCMC loop the function samples proposal values of the $\beta$ parameters from a MVN (the values for this function are set arbitrary to the values of beta0\_start, beta1\_start and sigma\_start. See Limitations) and it calls the logposterior() function to compute the posterior values at the proposal values and at the previous values of the cycle. Given the samplign of $log(\beta_1)$ adn the need for the exponent of this in the applied logistic model, the logposterior() function is called untill the value of computed log_ratio is different from NaN.
The function accepts the proposal folowing the rules of the MH algorithm. After the MCMC loop the diagnostics are collected for the chains: the acceptance rate (the tuning parameter of which has been set to 2.3 (to be checked)) and the Geweke z-values for the $\beta$ parameters. These diagnostics will be used at the end of the simulation to heck for the good convergence of the MCMC. 

The function takes in input the cohort dataset (cohort), the dataset witht he doses information (doses_info), the time by which the next patient is arrived (time_arrival), the simulation and the run values (i_simualtion and run) and the iterations (set to 10000) as well as the burn-in values (set to 1000). The function returns the posterior means of the $\beta$ parms and the diagnostic list. 

### MCMC_adaptive_EWOC()

The MCMC_adaptvie_EWOC() function is aimed to compute the MCMC posterior densities estimates for the $\beta_0$ and $\log(\beta_1)$ parameters. It does so applying an adaptive Metropolis-Hastings algorithm. 

Parameters: first $\lambda$ proposal is set to $log((2.38^2)/2)$ and the $\gamma$ is set to $\frac{1}{(1+i)^{\delta}}$, with $\delta = 0.01$. Accepetance rate = 0.24 in this way. However this approach has differences with the Nimble code (the non-adaptive version is far more simialr to the Imble one but it has lower than 24 ar). But most importantly the posterior distributons of the two parameters re too much equal one to the other. 
Type of adaptation: AM algorithm with global adaptive scaling (Ref: https://link.springer.com/article/10.1007/s11222-008-9110-y)
A small value is added to the updated version of Sigma to avoid being stuck anf following the idea in https://keichenseer.com/post/a-metropolis-algorithm-in-r-part-2-adaptive-proposals/

References: https://journals.sagepub.com/doi/pdf/10.1177/1536867X1401400309#cite.AT08@-11, https://biodatascience.github.io/statcomp/advmcmc/advmcmc.html#24_Adaptive_Metropolis_Algorithm, https://keichenseer.com/post/a-metropolis-algorithm-in-r-part-2-adaptive-proposals/


### MCMC_adaptive_EWOC()
The MCMC_adaptvie_EWOC() function is aimed to compute the MCMC posterior densities estimates for the $\beta_0$ and $\log(\beta_1)$ parameters. It does so applying an adaptive Metropolis-Hastings algorithm. 
Parameters: first $\lambda$ proposal is set to $log((2.38^2)/2)$ and the $\gamma$ is set to $\frac{1}{(1+i)^{\delta}}$, with $\delta = 0.01$. Accepetance rate = 0.24 in this way. Better similarity with Nimble ut most importantly the two distributons are not the same moved on the x-axis. 
Type of adaptive algorithm: Componentwise AM with componentwise adaptive scaling (Ref: https://link.springer.com/article/10.1007/s11222-008-9110-y)
A small value is added to the updated version of Sigma to avoid being stuck anf following the idea in https://keichenseer.com/post/a-metropolis-algorithm-in-r-part-2-adaptive-proposals/


### MCMC_nimble()

The MCMC_nimble() function is aimed to compute the MCMC posterior density estimates for the the $\beta_0$ and $\log(\beta_1)$ parameters. It does so applying an adaptive Metropolis-Hastings algorithm exploiting the packge Nimble. This is used just for comaprative purposes. Note that the reference_dose must be manually chnaged in the function specification. 

### decision()

The decision() function calls the MCMC() function to compute the updated \beta parameters. It then computes the probability of DLT at each of the doses and evaluates which is closest to the target of 0.3 (Use of the LRmean 'approach' by Neuen). The logistic model used is the following: $Logit(\pi_i) = \beta_0 + \beta_1d_i $, with $d_i$ the dose of subject i.

### evaluate_posterior()

The evaluate_posterior() is a function that is aimed to compute the posterior probability ranges as described by Neun. It takes in input a vector of posterior distributions and gives in output the probability that maximizes the targeted toxicity range while controlling for the EWOC for the excessive or unacceptable toxicity ranges. 
The function loops over all the computed a posteriori intervals and returns the probability index that maximizes it. It also checks if the potential probability is contolled for the EWOC (P(x >=0.35) > y, y = 0.25): if the selected one is not okay with the toxicity check then the loop is entered again to find the second to maximum targeted interval. The loop ends either with a valid dose or with a NaN.  
The function takes in input a limit dose that can resize the vector of possible doses to be evaluated (applied with the backfilling).

NOTE: the EWOC is defined to control the p(overtoxicity) to be less and equal to 25%.

### compute_parameters()

The compute_parameters() function is aimed to retrieve the parameters of a Beta distribution from the 25%, 50% and 75% quantile of it. It computes them starting from the distributions of the probabilities rv (Beta) an it involves the function rriskDistributions::get.beta.par(). 
It returns the dataframe with the Beta parameters for all of the evaluated probabilties.

### decision_EWOC()

The decision_EWOC is a function that is aimed to take the decision regarding the next dose. It calls the MCMC() and it uses the computed $\beta$ to compute the posteriori probabilities. It then uses the function evaluate_posterior() to retrieve the index of the most correct probability according to Neuen. It returns the choice (the next dose), the diagnostics from the MCMC and the data frame with the parametrs of the posterior Beta distributions (to be used as priors, info_doses).


### k_fold_skipping()

The k_fold_skipping() function checks the fold rise in the change of the dose: if the next dose is higher than k+1 times the current dose then we have to take the maximum affordable. 
Threshold for k-fold skipping has a default value of 2, meaning that the next dose can not be higher than 3 times the current dose.
2-fold rise means that A is 3 times B (2 folds more) (reference https://en.wikipedia.org/wiki/Fold_change) 

### hard_safety()

Check the hard_safety rule for all the doses used up to now in the cohort. Check with the constrain of the time_arrival of the pts that broke the backfilling cohort (you are interested in the data up to that arrival. This is for the backfilling pts that might not have finished the DLT window. Cohort are always ok by definiton so we can pass the time_arrival = 0 for them). 

The hard_safety() function checks whether at each dose the limit of maximum toxicity has been crossed. The limit is defined as the P(pdose_pDLT > \phi) \ge \psi, with pdose_pDLT the posterior probability at that dose. According to the paper (BF+BLRM) and to the BOIN setting, this is equivalent to P(pdpse_pDLT > 0.3) \ge 0.95. For this reason the same boundaries for eliminating a dose used in the BOIN setting are employed. The dose for which the inequality is true is reported as the limiting dose. 

### check_n_cap 

The check\_n\_cap() function checks if a dose should be closed for backfilling due to the n\_cap limit. It updates the doses\_info dataset if a dose should be closed (reference code = 2).
The function takes in input the cohort dataset and the doses info dataset and returns the updated doses\_info dataset

### transform_data()

The transform\_data() function is only aimed to properly set the class of the cohort data frame columns. It takes in input the original cohort dataset and it returns the transformed one. 


### transform_data_doses()

The transform\_data() function is only aimed to properly set the class of the cohort data frame columns. It takes in input the original cohort dataset and it returns the transformed one. 

### compute_BFBLRM()

The compute\_BFBOIN is a function devoted to performing the allocation procedure following the BF+BLRM logic. The functions takes in input the cohort dataset (cohort), the dataset of doses info (doses_info), the dataset with the boundaries for the decision (new\_reference), the limiting values for the total patients, the total patients for backfill and the total patients for cohort (n\_max, n\_cap and n\_stop), the time of followup (DLT\_time), the rate of accrual during the follow-up period (lambda), the target toxicity rate (target, default = 0.3) and the simulation parameter (i_simulation, default = 0). It returns the cohort dataset with all the patients allocated and the updated dataset with the diagnostics matrix (to be updated with the safety and sufficeint info values). 

The function works with a while loop that is entered until the n\_max is reached, or the minimal dose shows overtoxicity or the dufficient info is reached (see Appendix 1). The parameters set prior to the while loop are: 
- run: the counter for the runs (a run corresponds to the beginning of a new cohort)
- pts: the counter for the number of patients in the trial
- n\_early\_stop: counts the number of times a trial is stopped due to overtoxicity (simulation)
- n\_stop\_reached: counts how many times a trials stopped due to the sufficient info rule 
- last\_closed_dose: a variable indicating the last closed dose for the backfilling (set to NA at the beginning) 
- limit\_time: variable for the cohort() function, set to 0 at the beginning
- current\_dose: the current dose used in the trial. Set to the minimum dose at the beginning of the trial
- last\_arrival: set to 0 see backfill\_patients() and cohort\_patient()
- last\_lag: set to 0 see backfill\_patients() and cohort\_patient()
- doses: the vector of doses from the doses\_info dataset. Used just for ease
- diagnostics: a matrix of 3 columns to store the diagnostics of the MCMC 


The while loop has two mutually exclusive sections: one for the case in which no backfill doses has been opened yet (state vector is a NA vector) and one for when the backfilling has been opened. The opening of backfilling depends on whether at least one response has been seen and on whether an escalation from the first dose has been made. To enter one of the two sections, the boolean variable 'backfill' is used. The default value of it is T (allowing for the possibility backfilling) and this value is set to F (not allowing the possiblity of backfilling) if the next dose is the minimum one or the current dose is the maximum one. If backfill = T, then the probability of response at the current dose is simulated and if at least one response is observed, the 'State' of the doses is updated (0 to all and 1 to the current observed dose. It might be that the currently open dose is not the lowest but the opening of lower doses is forbidden as no response has been seen there). The current dose is then deemed safe and open to backfilling. The doses higher than the current one are potentially open in term of probability of response. If the current dose is not open since of no response then it is set to -1.

The first passage in each of the sections is to compute the current cohort patients values. In order to do this the function cohort\_patients() is used.
After the function call, the total number of DLT at the current dose and the total number of patients at the current dose are computed. There is no need to check whether all the patients are finished, since we are just following this cohort and we assume that after the cohort\_patient() function, all patients have been followed-up. Since the backfilling is not open yet there is not need to check for those patients. The computed values are used to take the decision that is computed with the decision() function. The result is stored for check purposes. The decision is then evaluated in terms of the K-Fold skipping rule, the Hard safety rule and the Sufficient Information rule (see Appendix 1). IF a limit dose is found all the doses above (including the limiting dose) are closed. If the limiting dose is the minim dose then the trail is stopped for safety. IF the next dose is the same as the current dose then the Sufficient Information rule is checked and eventually the trial is stopped.

If there is a possibility for backfilling (backfill = T), the responses at the current dose level are computed. Up to now there are two strategies to do so: 
- 1) use of a binomial with n = the cohortsize under investigation and p = the probability of response at the selected dose. If the binomial is $\ge 1$ then the current dose is open for backfilling, with the state set to 1 and all the other states set to 0
- 2) use of the Weibull distribution, with DLT\_time the quantile at the probability of response value. For each of the new patients the value of the time for the response is computed an if it follows within the patients' limit\_time (maximum follow-up time, that corresponds to the maximum of the inidividual limit times) then the response is 1 otherwise it is 0 (responses after the limit\_time are not counted). If the number of responses at the current dose is $\ge 1$ then the state of the current dose is updated to 1 and all the others to 0. The limit time is computed only across the current run patients.

Suppose the backfill has been opened.The cohort patients are assigned thanks to the cohort\_patient() function and the maximum time of follow-up across the actual cohort patients is computed (limit\_time). The arrival of other patients is simulated: all those patients that arrive between the time of arrival of the last cohort patient (previous\_time as returned by the cohort\_patient() function) and the limit\_time must be put in a backfill arm. The process is handled by means of the backfill\_patients() function.

After having described the arrival of all the patients, all the new patients are assigned to a specific dose level. Starting from the earliest patient, the maximum\_open() function is called to find the appropiate dose level.If there is at least one dose open then it is checked if this should be closed due to new information (patients ending the limit\_time by the arrival of the current subject). This process is aimed to eventually change the current maximal dose and it is doen with the HArd safety rule (see Appendix 1 and Limitations).
In the case where no dose is found open by the arrival of the patient (maximum\_open() returns NA) and the current backfill dose is NA then the patient is excluded from the trial (no info update) and the loop moves to the next patient.

In the case the patient is included, his/her infos are then updated in the appropriate manner. Also, the check\_n\_cap() function is called to check whether after the insertion of the new patient the dose is still potentially open for backfilling (state different from 2).

After having assigned all the patients in the backfill cohort, a filtering procedure is undergone to remove those patients that were not assigned to either dose as by their arrival no open dose was found.
After all the assignements, the merged\_decision() is invoked and the next dose for the cohort is found. The same rules explained before are applied to evaluate the decision for the next dose, except for the backfill variable that is not here used.\\
After the dose assignment, the current dose is checked to eventually change its backfill status from -1 to 1 (check for activity, with the same procedure previsously explained).

The backfill doses are also checked to close all those above the next dose and to evetually open all those below the next dose (unless a code for enough patients or not activity is present. See Appendix 3). 

The current dose (for the next cycle) is updated with the next dose found from either the branches of the main while loop. 

Note that in case the cohort has been assigned to a previously closed dose for backfill (due to no response), the evaaluation for the opening of this is done again and eventually the state is altered to 1 (the dose won't be used for backfilling if the cohort is at the same dose and this is due to the checkings of <= current_dose).

THe function returns the entire patient set up dataframe (cohort), the counter (binary) for the early stops for safety and the counter (binary) for early stop for sufficient information.


### compute_BFBLRM_EWOC()

This function is analogous to compute_BFBLRM() one but it applies the rules of Neuens. for the finding of the optimal dose and it applies a different metholodogy to find deal with the backfilling.

### BFBLRM()

The function BFBRLM() takes in input the doses, the true probabilities of DLT and response, the cohortsize, the maximum number of patients, the total patients for backfill and the total patients for cohort (n\_max, n\_cap and n\_stop), the DLT time (DLT\_time), the rate of accrual during the DLT time (lambda) and the target DLT rate (target). It creates the dataframe doses\_info with all the information regarding the dose levels (dose, probability of DLT and response, the shape and scale parameters of the Weibull distributions and the state). It calls the compute\_BFBLRM() to compute the allocation of the patients. It returns then the object results.

## Appendix

### 1 BLRM + BF RULES (to be finisehd, chosen rules of BF+BLRM paper)
- Hard Safety
- K-fold skipping 
- Sufficient information
 
### 2: Time values 

- Lag and lag_arrival: distance measure between the last patient arrival and the current patient arrival (within the same cohort call or backfill call. Otherwise it is the difference between the limit time and the next patient arrival)
- Time_arrival and time_arrival: a point in the trial timeline that indicates when the current patient has arrived. Computed as the sum of the previous patient arrival point and the lag of the current patient
- Time_DLT and time_DLT: a distance measure between the arrival point of the current patient and his/her DLT experience
- Limit_time: a point in the trial timeline indicating when the current patient is stopped to be followed-up. Computed as the minimum between the patient's time of arrival summed to the DLT_time and the patient's time of arrival and the time of the DLT experience. In the cohort_patient() function it is used as the previous_time point from which the next cohort patients' arrival times are generated. in the backfill_patients() function it is used as the maximum time point behind which the backfill can be used
- DLT_time: follow-up time equal for all the patients. Distance measure between the arrival of a patient and the forced end of the follow-up
- last_arrival_time: a point in the timeline that indicates the arrival of the first patient after the backfilling period
- last_lag: a distance indicating the lag between the last arrival time and the first no-backfill patient's arrival time
- Previous_time: wihitn the functions of cohort_patient() and backfill_patients() represents the last time point in the time line. For the cohort_patient() function: the input value is he previous cohort limit time (regardless of a backfill or not) such that all the new cohort patients arrives after the maximum follow-up time; the output value is the last patient's arrival time. In the backfill_patients() function: the input value is the last cohort patient's arrival time so that the backfill patients arrive between the last patient's arrival point and the maximum follow-up of the cohort; the output value is the last backfill patient's arrival time (not necessary to be returned)
- Sum_times: a point in the timeline that indicates when a patient has a DLT (sum of arrival time and time_DLT)

### 3: Codes for backfill
- 0: closed
- 1: open 
- 2: closed for n\_cap
- -1: closed for not activity

## Limitations

-  In the decision() function we report the selected probability as the one that minimizes the absolute distance with the target. In the case for which two are minimizing then the lowest one is taken. Can be a conservative approach
-  In the lowest-Dose-deemed-unsafe part I use the same limit as in the BOIN (0.95) whereas in the BLRM paper this limit is set to 0.80 (more conservative).
-  Backfill closed with the overtoxicity boundary (different from BOIN)
-  The starting values in the MCMC are arbitrarly chosen as for the few observations a var-cov matrix estimated in any way would be extremelly large and so then the loop for the NA of the log_ratio would be never exit

# What to do 9/4/24: 
- finish commenting
- discuss problems
- check code Neuen
- finish with implementation of the simulation settigns 
