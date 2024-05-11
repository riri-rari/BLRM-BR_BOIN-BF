This folder contains the results of the simulations done with the BFBOIN.R and BFBOIN_at.R. 
The scenario run are 

### Scenario 0
true_pDLT <- c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55, 0.6)

# Scenario 1
true_pDLT <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
true_presp <- c(0.3, 0.4, 0.45, 0.5, 0.55, 0.6)

# Scenario 2
true_pDLT <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
true_presp <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# Scenario 3
true_pDLT <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
true_presp <- c(0.1, 0.2, 0.3, 0.45, 0.58, 0.67)

# Scenario 4 
true_pDLT <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
true_presp <- c(0.05, 0.1, 0.15, 0.3, 0.45, 0.55)

# Scenario 5
true_pDLT <- c(0.05, 0.1, 0.15, 0.20, 0.3, 0.4)
true_presp <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)

#scenario 6
true_pDLT <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.3)
true_presp <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)

# Scenario 7 (sharp change with constant prob of response)
true_pDLT <- c(0.05, 0.05, 0.05, 0.8, 0.8, 0.8)
true_presp <- c(0.3, 0.3, 0.3, 0.36, 0.36, 0.36)

# Scenario 7 (sharp change with constant prob of response)
true_pDLT <- c(0.05, 0.05, 0.05, 0.8, 0.8, 0.8)
true_presp <- c(0.3, 0.3, 0.3, 0.36, 0.36, 0.36)

# Scenario 8 (non linear relationship around the MTD with precise target)
true_pDLT <- c(0.15, 0.2, 0.25, 0.3, 0.45, 0.6)
true_presp <- c(0.15, 0.3, 0.35, 0.36, 0.38, 0.40)

# Scenario 9 (non linear relationship around the MTD with precise target)
true_pDLT <- c(0.05, 0.15, 0.3, 0.35, 0.4, 0.45)
true_presp <- c(0.10, 0.2, 0.3, 0.35, 0.4, 0.45)

# Scenario 10 (equal to scenario 2 but with different p_responses)
true_pDLT <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
true_presp <- c(0.1, 0.15, 0.2, 0.25, 0.30, 0.35)

# Scenario 11 (equal to scenario 5 but with different probabilities)
true_pDLT <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)
true_presp <- c(0.05, 0.1, 0.12, 0.17, 0.2, 0.25)

Constant values: 
doses <- seq(1, 6, 1)
cohortsize <- 3
n_cap <- 12
n_stop <- 9
n_max <- 32
DLT_time <- 1 #1 month
lambda <- 6 
