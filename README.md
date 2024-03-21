# Introduction 

This code is aimed to implement the BF+BOIN algorithm used in Phase I oncology CT. 

## Backfill + BOIN 

### Baseline functions 
#### weibull_parms()

weibull_parms is a function that takes in input the probability of an event (as DLT at a specific dose, ie pDLT) and the time within which the event might happen (DLT_time, that is the time of patient follow-up). The function assumes that the DLT_Time is the pDLT quantile of a Weibull function and that $0.5*DLT_time$ is the $0.5*pDLT$ quantile of a Weibull distribution. From the quantiles, the function computes the values of the Shape and Scale parameters (ref: https://www.johndcook.com/quantiles_parameters.pdf).

