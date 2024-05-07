# Introduction 

This code is aimed to implement the BF+BOIN algorithm.

## Backfill + BOIN 

### Baseline functions 
#### weibull_parms()

weibull_parms is a function that takes in input the probability of an event (as DLT at a specific dose, ie pDLT) and the time within which the event might happen (DLT_time, that is the time of patient follow-up). The function assumes that the DLT_Time is the pDLT quantile of a Weibull function and that $0.5*DLT_time$ is the $0.5*pDLT$ quantile of a Weibull distribution. From the quantiles, the function computes the values of the Shape and Scale parameters (ref: https://www.johndcook.com/quantiles_parameters.pdf).

#### cohort_patient()

The cohort_patient function is devoted to the computation of the values for the cohort patients. It takes in input  the cohort dataset (cohort), the info about the doses (doses\_info), the current dose value (current\_dose), the current run (run), the current patient count (pts), the current previous time (previous_time), the last arrival time (last\_arrival\_time, default = 0 otherwise eqaul to the time of the person who stopped the backfilling. See backfill\_patients()) and the last lag (last\_lag, default = 0 otherwise set to the last lag of the person who stopped the backfilling. See backfill\_patients()). The function returns the cohort dataset updated with the patients info, the updated patients count and the updated previous time.

The function retrieves the parameters of the Weibull distribution at the current dose from the doses\_info dataset. Then, it cycles over a cohortsize of patients: if the patient is the first one of the trial, then his/her time\_arrival and lag\_arrival are not computed (set to 0) and the previous time is kept to the original value (it should be passed as 0); if the patient is not the first one, then his/her lag of arrival is simulated from the Exp(lambda), the time\_arrival is computed as the sum of the previous\_time (that represents the time of arrival of the previous patient. The input value of previous\_time is the limit\_time of the previous cohort) summed to the lag\_arrival and the previous\_time is then updated to the current patient time\_arrival. If a value for the last\_arrival\_time is given and the call of the function is 1 (so it is the first patient at the current cohort), the time\_arrival of the patient is set equal to the given arrival time and the lag is set equal to the given lag time (this is thought to reusing the patient that has broken the bakfilling process. See backfill\_patients()).

Regardless of the patient's number, the time for the DLT is computed (time\_DLT), the DLT is computed and the limit_time is computed (this represents the point in time when the patient is stopped to be followed-up).

The computed information are returned as a list. A summary of the time values is in Appendix 2

#### decision()

The decision function is devoted to consult the table of decisions (reference\_general. Computed with the BOIN::get.boundary() function or the get\_boundaries() function under the same priors scenario) at the current number of patients. The function then outputs the decision (change) that can take four values: 1 = 'ESCALATE', 2 = 'STAY', 3 = 'DESCALATE', 4 = 'ELIMINATE'. The function takes care of the review from FDA https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.fda.gov/media/155364/download&ved=2ahUKEwjDpImqjfuFAxVLSfEDHaEXBFgQFnoECBMQAQ&usg=AOvVaw2gMeg-rLC1ZMSKTZqx1dBa 

#### backfill_patients()

The backfill\_patients() function is aimed at the simulation of all those patients that should enter a backfill arm as they arrive between the last cohort patient and the limit time. It returns an updated version of the cohort dataset with the backfillable patients (info of their dose and pDLT-related characteristics are set to NA) and the time that has make the while loop ceasing (the over-limit-time-patient's arrival time). 

Time arrival constrain is less-equal to the limit_time. 
Patients are generated up to the n\_max limit if all the patients arrives within the limit time.

THe function resturns the updated cohort, the updated number of patients, the updated previous time (last time of arrival of a backfilled patient), the last time arrival (last\_arrival) and the last lag (last\_lag) corresponsing to the arrival time and the lag time of the patient that interrupted the backfill.

#### open_close()

The open_close() function is used to decide whether a specific dose should be open or closed for backfilling. This is done in in three steps: 
- a) collection of the total patients and DLT events at the same dose level. Only patients with (time_arrival + DLT_time) $\le$ current backfilled patient's arrival time are included (the only info that we have at this time point)
- b) collection of the total patients and DLT events at the dose level immediately higher than the current. Only patients with (time_arrival + DLT_time) $\le$ current backfilled patient's arrival time are included (the only info that we have at this time point)
- c) computation of the decision (decision() function) for a and a + b datasets and evaluation of the results. If the dose is closed for backfilling, all the doses greater than this are closed too (see maximum_open()). 
- 
If a dose should be closed then a 0 is returned, while if it should be open, function returns 1. In the case for which no patietns are available to have info by the arrival of the new patient, then the dose is open (return 1)

The function takes in input the cohort dataset, the time of the current patient of backfll (time\_pts\_backfill), the dose that is currently under investigation (dose), the reference dataset for the boundaries (new\_reference. See BFBOIN()) and the dataset of info about the doses (doses\_info).

The function returns ether 1 or 0, as specified above. 

#### maximum_open()

The maximum_open() function serves to find the maximum dose open for the backfill. The function loops over the doses (each time it takes the highest dose open for backfill and lower than the current cohort dose) until a dose can be kept open (change\_backfill switches from F to T). The function returns the found backfill dose and the updated doses\_info dataset.

The function is deemed to check if the previously highest dose open for backfill is still okay and to eventually close it. Since the function loops until a change of dose is not found, all those doses that were open for backfill prior to the arrival of the current patient and that by the arrival of the patient should be closed, are closed. The maximum dose open is then found. 
This takes into consideration all those patients who has experienced a DLT between the previous patient and the current patient arrival time, changing the evaluation of the doses for backfilling. In the case for which no more open doses are left for backfilling (current\_dose\_backfilling is an empty object), the while loop is broken and the function returns a NA value. In the case that the current dose is available, the state of it is checked by means of the open\_close() function. 
If this function returns 1 then the dose is left open, while if it returns 0 the dose is closed as well as all the doses that are higher than the current one. 

In the case that a dose is closed, then the cycle is entered again. If the dose is, instead, the right one (open\_close returns 1), the cycle is exit (change = T) and the current dose is returned. 

The function takes in input the cohort dataset, the time of the current patient of backfll (time\_pts\_backfill), the dose that is currently under investigation (current\_dose), the dataset of info about the doses (doses\_info) and the reference dataset for the boundaries (new\_reference. See BFBOIN()) 

#### check_n_cap 

The check\_n\_cap() function checks if a dose should be closed for backfilling due to the n\_cap limit. It updates the doses\_info dataset if a dose should be closed (reference code = 2).

The function takes in input the cohort dataset and the doses info dataset and returns the updated doses\_info dataset.

#### conflict()

The conflict() function is aimed to compute the maximum dose of the backfill for which there is a conflict between the backfill decision and the current cohort dose decision. To do this, the function cycles over all doses values that have been assigned during the backfilling process and it registers the maximum dose for whcih the conflict table is eqaul to 0 (see Appendix for codes). The function then returns the maximum dose found and the decision that was computed at this dose. In case of no doses assigned during the backfilling or no conflict the maximum dose is set to NA and the change is set to NA (in the merged\_decision() fucntion, this will be read as no need to look at the backfilling decision but just at the cohort decision). 

The function takes in input the cohort dataset, the current cohort dose (dose\_cohort), the new\_reference table, the list of actual run doses of the backfill (dose\_back\_run) and the value of the arrival of the next patient (last\_arrival). 

#### merged_decision()

The merged\_decision() function is aimed to compute the decision between the backfill and the cohort datasets. It returns the safe dose found. It calls the decision() function inside and, since it is provided with the conflict table (ref. conflict table of the paper 'Backfilling Patients in Phase I Dose Escalation Trials Using Bayesian Optimal Interval Design (BOIN)'), it can handle the cases in which the two decisions are not agreeing. In the case of conflict, the cumulative ratio (q\_k) is computed from the value of b (backfill dose) to the value of k (this is the cohort dose at the beginning). Following the rules of the paper, if the current cumulative q\_k is above the lambda\_2 (threshold for the safety) then the highest dose j with q\_j below the threshold (with $b \le j \le c-1$) should be taken as the safe dose. To achieve this, the function computes the decision for the cumulative data from the highest dose between b and c and stops when a safe dose is found (if the dose is the maximum dose and the decision is to stay at the current case. This scenario should be never reached as the current dose is at most c-1). In the case that no safe dose exists or it is found between b and c-1 then the function assigns the b-1 dose as the safe one (in the case the backfill dose is already the minimum dose then the safe dose is left to be the backfill dose).

In the case of no conflict or in case no backfill dose info is available (dose_backfill$max_dose = NA), the decision from the cohort is kept and the safe dose is assigned according to the rules of BOIN. 

The function returns both the safe dose and the decision, so that the evaluation of a closure for toxicity or for maximum number of people (+ 'STAY') can be evaluated.

Note that the patients from the backfill dose are selected according to the last\_arrival value: those are all the patients that can contribute to the decision by the arrival of the new cohort patient (same value that made the backfill\_patient() loop breaks). This value can be changed to the limit\_time of the cohort arm (implying that the backfill patients whose DLT is visible by the end of the cohort follow-up are used. This possibly leads to a loss of info from all the patients that completes between the limit\_time and the arrival of the new patient).

### transform_data()

The transform\_data() function is only aimed to properly set the class of the cohort data frame columns. It takes in input the original cohort dataset and it returns the transformed one. 

### transform_data_doses()

The transform\_data() function is only aimed to properly set the class of the cohort data frame columns. It takes in input the original cohort dataset and it returns the transformed one. 

### compute_BFBOIN()

The compute\_BFBOIN is a function devoted to performing the allocation procedure following the BF+BOIN logic. The functions takes in input the cohort dataset (cohort), the dataset of doses info (doses_info), the dataset with the boundaries for the decision (new\_reference), the limiting values for the total patients, the total patients for backfill and the total patients for cohort (n\_max, n\_cap and n\_stop), the time of followup (DLT\_time), the rate of accrual during the follow-up period (lambda) and the simulation parameter (i_simulation). It returns the cohort dataset with all the patients allocated, the safety and the sufficient info values. 
The function works with a while loop that is entered until the n\_max is reached, or the minimal dose shows overtoxicity or the n\_stop + 'STAY' decision is reached. The parameters set prior to the while loop are: 
- run: the counter for the runs (a run corresponds to the beginning of a new cohort)
- pts: the counter for the number of patients in the trial
- n\_early\_stop: counts the number of times a trial is stopped due to overtoxicity (simulation)
- n\_stop\_reached: counts how many times a trials stopped due to the c_stop + 'STAY' (simulation)
- last\_closed_dose: a variable indicating the last closed dose for the backfilling (set to NA at the beginning) 
- limit\_time: for cohort() and set to 0 at the beginning 
- current\_dose: the current dose used in the trial. Set to the minimum dose at the beginning of the trial
- last\_arrival: set to 0 see backfill\_patients() and cohort\_patient()
- last\_lag: set to 0 see backfill\_patients() and cohort\_patient()
- doses: the vector of doses from the doses\_info dataset. Used just for ease
- previous\_time: for the cohort() function and set to 0 

  
The while loop has two mutually exclusive sections: one for the case in which no backfill doses has been opened yet (state vector is a NA vector) and one for when the backfilling has been opened. The opening of backfilling depends on whether at least one response has been seen and on whether an escalation from the first dose has been made. To enter one of the two sections, the boolean variable 'backfill' is used. The default value of it is T (allowing for the possibility backfilling) and this value is set to F (not allowing the possiblity of backfilling) is a 'STAY' decision is made at the lowest dose, if the current dose receives a 'DESCALATION' or 'ELIMINATE' decision or if the curent dose is the maximum dose (Appendix 1). If backfill = T, then the probability of response at the current dose is simulated and if at least one response is observed, the 'State' of the doses is updated (0 to all and 1 to the current observed dose. It might be that the currently open dose is not the lowest but the opening of lower doses is forbidden as no response has been seen there). The current dose is then deemed safe and open to backfilling. The doses higher than the current one are potentially open in term of probability of response. If the current dose is not open since of no response then it is set to -1.

The first passage in each of the sections is to compute the current cohort patients values. In order to do this the function cohort\_patients() is used.
After the function call, the total number of DLT at the current dose and the total number of patients at the current dose are computed (points a and b). There is no need to check whether all the patients are finished, since we are just following this cohort and we assume that after the cohort\_patient() function, all patients have been followed-up. Since the backfilling is not open yet there is not need to check for those patients. The computed values are used to take the decision that is computed with the decision() function (point c). The result is stored for check purposes. The decision is then evaluated in relation with the current dose (point d): 
- the current dose is the minimum one and the decision is the 'ELIMINATE', the trial is stopped (backfill = F)
- the current dose is not the minimum one and the decision is to 'ELIMINATE', then all the doses $\ge$ current dose are eliminated and the next dose is the one below the current dose (backfill = F)
- the current dose is not the highest dose and the decision is to 'ESCALATE', then the next dose is the next dose (backfill = T)
- the current dose is the maximum dose and the decision is to 'ESCLATE', then the next dose is the current dose (backfill = T) and the check for n\_stop() is done
- the current dose is not the minimum dose and the decision is to 'DESCALATE', then the next dose is the lower dose (backfill = F)
- the current dose is the minimum dose and the decision is to 'DESCALATE', then the next dose is the current dose (backfill = F) and the check for n\_stop is done
- the current dose is not the minimum dose and the decision is to 'STAY', then the next dose is the current one (backfill = T)
- the current dose is the minimum one and the decision is to 'STAY', then the next dose is the current one (backfill = F)
- the decision is to 'STAY' so the current dose is checked for the n\_stop. 

Eventually the n_stop_reached counter is updated and the while loop is stopped as well as the trial
In the case the current dose is the maximum one, then the backfill is set to F (since dose escalation done without backfilling open so no dose prior to this had responses/activities and the actual dose is the maximum so it can not be used to have backfilling). 

If there is a possibility for backfilling (backfill = T), the responses at the current dose level are computed. Up to now there are two strategies to do so: 
- 1) use of a binomial with n = the cohortsize under investigation and p = the probability of response at the selected dose. If the binomial is $\ge 1$ then the current dose is open for backfilling, with the state set to 1 and all the other states set to 0
- 2) use of the Weibull distribution, with DLT\_time the quantile at the probability of response value. For each of the new patients the value of the time for the response is computed an if it follows within the patients' limit\_time (maximum follow-up time, that corresponds to the maximum of the inidividual limit times) then the response is 1 otherwise it is 0 (responses after the limit\_time are not counted). If the number of responses at the current dose is $\ge 1$ then the state of the current dose is updated to 1 and all the others to 0. The limit time is computed only across the current run patients.
  3) 
Suppose the backfill has been opened.The cohort patients are assigned thanks to the cohort\_patient() function and the maximum time of follow-up across the actual cohort patients is computed (limit\_time). The arrival of other patients is simulated: all those patients that arrive between the time of arrival of the last cohort patient (previous\_time as returned by the cohort\_patient() function) and the limit\_time must be put in a backfill arm. The process is handled by means of the backfill\_patients() function. 

After having described the arrival of all the patients, all the new patients are assigned to a specific dose level. Starting from the earliest patient, the maximum\_open() function is called to find the appropiate dose level.If there is at least one dose open then it is checked if there are openable doses between the current backfill dose and the cohort dose (excluding them). This process is aimed to eventually change the current maximal dose open to another one: if by the arrival of the new patient a previously closed dose (0 state) can be re-opened, then the patient should be assigned to this dose (since each time the patient should be assigned to the highest dose open for backfilling).
To do so, a vector with all the doses between the current backfill dose and the current cohort dose is created (higher\_doses; only doses with state different from 2 are taken) and if it contains at least one element a while loop is entered too. The while loop calls the open\_close() function for each of the doses in higher\_doses: if the current dose can not be re-opened then all the higher doses with respect to this can not be opened too and the loop is exit; if the dose can be open, the loop is entered again to test the next dose. In the case where no dose is found open by the arrival of the patient (maximum\_open() returns NA), the higher\_dose vector is started from the minimum dose with response (included) to the current cohort dose (excluded). At the end, the current backfill dose is update with the appropriate dose and if the dose is NA then the patient is excluded from the trial (no info update) and the loop moves to the next patient.

In the case the patient is included, his/her infos are then updated in the appropriate manner. Also, the check\_n\_cap() function is called to check whether after the insertion of the new patient the dose is still potentially open for backfilling (state different from 2)

After having assigned all the patients in the backfill cohort, a filtering procedure is undergone to remove those patients that were not assigned to either dose as by their arrival no open dose was found (not even re-openable).

After all the assignements, the merged\_decision() is invoked and the next dose for the cohort is found. The same rules explained before are applied to evaluate the decision for the next dose, except for the backfill variable that is not here used.\\
The current dose (for the next cycle) is updated with the next dose found from either the branches of the main while loop. 

In the case of a different choice than 4 for the next dose, the current dose is checked for eventually open the backfill in the case of closing for previous inactivity. If the dose has already reached the limit number of patients then the code is not changed.

Note that in case the cohort has been assigned to a previously closed dose for backfill (due to no response), the evaaluation for the opening of this is done again and eventually the state is altered to 1 (the dose won't be used for backfilling if the cohort is at the same dose and this is due to the checkings of <= current_dose). 

THe function returns the entire patient set up dataframe (cohort), the counter (binary) for the early stops for safety and the counter (binary) for early stop for sufficient information.

### BFBOIN()

The function BFBOIN() takes in input the doses, the true probabilities of DLT and response, the cohortsize, the maximum number of patients, the total patients for backfill and the total patients for cohort (n\_max, n\_cap and n\_stop), the DLT time (DLT\_time), the rate of accrual during the DLT time (lambda) and the target DLT rate (target). It creates the dataframe doses\_info with all the information regarding the dose levels (dose, probability of DLT and response, the shape and scale parameters of the Weibull distributions and the state). It calls the get.boundary() function from the BOIN library (can be substituted with the other one) and then it calls the compute\_BFBOIN() to compute the allocation of the patients. It returns then the object results. 

## Implementation 

Description: 
- doses: an increasing-ordered vector of doses to be tried in the trial 
- true_pDLT: a vector of true probabilities of toxicity at the dose levels. The order reflects the order in 'doses'
- true_presp: a vector of true probabilities of response at the dose levels. The order reflects the order in 'doses'
- cohortsize: the size of the group of patients to be added at each cohort
- n_cap: the number of patients that stops the trial if a 'STAY' decision follows at the same dose
- n_stop: the number of patients above which backfilling at the same dose is forbidden
- n_max: the maximum number of patients in the study
- DLT_time: the follow-up time (for now specified in months)
- lambda: the rate of accrual people in a DLT_time
- doses_info: a data frame with the information regarding the probability of DLT, the Shape and the Scale of the Weibull(DLT), the probability of response, the Shape and the Scale of the Weibull(Response) and the state for each dose level
- cohort: the dataframe for the trial patients. It includes the Run at which the patient has been included (a run is a cycle in the while loop. A run changes with a cohort change), the Patient number (increasing from 1 to n_max), the group of belonging ('C' for cohort, 'B' for backfill), the dose at which the patient is exposed, the time at which the patient has arrived ('Time_arrival'. The time starts with the first patient that has 'Time_arrival' set to 0), the time of the DLT ('Time_DLT', it indicates when the patient experiences the DLT. It is a relative time with respect to Time_arrival), the indicator for having a DLT or not (0 is Time_DLT is above DLT_time and 1 if it is below), the distance in time with respect to the previous patient ('Lag', coming from Exp(lambda)) and the maximum time a patient is followed-up (computed as the minimum between (Time_arrival + DLT_time) and (Time_arrival + Time_DLT). This reflects the idea that if a patient experiences a DLT before the ending of the DLT_time then the follow-up time is interrupted). The summary of the time values is in Appendix 2

### Appendix 

#### 1: Cases for the backfill variable 

If the current dose is the minimum one: 
- 'ELIMINATE' = terminate the trial, no problem with backfilling
- 'DESCLATE' = 'STAY' = no backfilling
- 'ESCALATE' = ok possibility of backfilling

If the current dose is not the minimum one: 
- 'ELIMINATE' = no backfilling at that dose
- 'DESCALATE' = no backfilling at that dose
- 'STAY' = ok possibility of backfilling (further chacked in later passages)
- 'ESCALATE' = ok possibility of backfilling 

If the current dose is the maximum one: no need to check for backfilling 

#### 2: Time values 

- Lag and lag_arrival: distance measure between the last patient arrival and the current patient arrival (within the same cohort call or backfill call. Otherwise it is the difference between the limit time and the next patient arrival)
- Time_arrival and time_arrival: a point in the trial timeline that indicates when the current patient has arrived. Computed as the sum of the previous patient arrival point and the lag of the current patient
- Time_DLT and time_DLT: a distance measure between the arrival point of the current patient and his/her DLT experience
- Limit_time: a point in the trial timeline indicating when the current patient is stopped to be followed-up. Computed as the minimum between the patient's time of arrival summed to the DLT_time and the patient's time of arrival and the time of the DLT experience. In the cohort_patient() function it is used as the previous_time point from which the next cohort patients' arrival times are generated. in the backfill_patients() function it is used as the maximum time point behind which the backfill can be used
- DLT_time: follow-up time equal for all the patients. Distance measure between the arrival of a patient and the forced end of the follow-up
- last_arrival_time: a point in the timeline that indicates the arrival of the first patient after the backfilling period
- last_lag: a distance indicating the lag between the last arrival time and the first no-backfill patient's arrival time
- Previous_time: wihitn the functions of cohort_patient() and backfill_patients() represents the last time point in the time line. For the cohort_patient() function: the input value is he previous cohort limit time (regardless of a backfill or not) such that all the new cohort patients arrives after the maximum follow-up time; the output value is the last patient's arrival time. In the backfill_patients() function: the input value is the last cohort patient's arrival time so that the backfill patients arrive between the last patient's arrival point and the maximum follow-up of the cohort; the output value is the last backfill patient's arrival time (not necessary to be returned)
- Sum_times: a point in the timeline that indicates when a patient has a DLT (sum of arrival time and time_DLT)

#### 3: Values for the decision 

- 4: Elimination 
- 3: De-escalate
- 2: Stay
- 1: Escalate

#### 4: Codes for conflict 

- 0: conflict 
- 1: no conflict

#### 5: Codes for backfill 
- 0: closed
- 1: open 
- 2: closed for n\_cap
- -1: closed for not activity

#### 6: Limitations 

- If BY the arrival of a patient no dose is open for backfill and no dose is openable for backfill then the patient is rejected 
- In case of conflict and the backfill dose is the minimum dose and there is no safe dose between backfill and cohort, then the backfill is assigned as safe dose 
- If there is no a cohortsize number of spots available for the next cohort, the algorithm returns
- Once i have the dose escalation then i need to open the backfill wt dose-1 without caring for safety (makes sense as escalation done + No problem with re-opening since for a re-opening you have for sure patients  ). But then i do not care of how the actual dose is in terms of safety: should i check if the dose is still good? even if i do not have an info from same dose or upper dose? Bcs i do have patients --> you have checked with the the open_close and in case of no pts available then you do not alter the current level (equal to say that you do not have data to update/alter the info)
- The DLT_time is used as the time for assignment (reference to BLRM and TITE-BLRM)
- Titration might lead to a lot of underdosing ppl for the dose escalation --> I think it is the same in the paper (+ they consider all the ppl for underdosing of overdosing)

#### 6.bis: Scenari

- multiple doses as those from the BOIN paper
- change of the target? reasonable? 
- change of accrual rate? Change of DLT window? 
- change of the safety boundary (for the elimination)?





