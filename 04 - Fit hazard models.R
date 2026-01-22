# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04 - Fit hazard models
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 25 Nov 2025
# Date last modified: 22 Jan 2026 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(nimble)          # modeling with nimble
library(nimbleHMC)       # HMC in nimble
library(coda)            # model outputs
library(mgcv)            # construct basis functions for cyclic splines

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_forModel.csv")

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# create the basis matrix (cyclic spline)
basis <- cSplineDes(fates$week,
                    knots = knot.list,
                    ord = 4)                # cubic

#_______________________________________________________________________________________________
# 4. Model setup ----
#_______________________________________________________________________________________________
# 4a. Build data lists ----
#_______________________________________________________________________________________________

# constants
constants = list(
  
  N = nrow(fates),              
  nsf = 4,               # sex/forest type combinations
  nbasis = ncol(basis),
  
  # constant indices
  s_f = fates$sex_forest
  
)

# data
fates.1 <- list(
  
  # responses, split by scenario
  y_mort_1 = fates$y.mort.scen1,
  y_mort_2 = fates$y.mort.scen2,
  y_mort_3 = fates$y.mort.scen3,
  
  # covariates
  bci_1 = (fates$BCI.1 - mean(fates$BCI.1)) / sd(fates$BCI.1),   # standardized (center + scale)
  post1 = fates$post1,
  post2 = fates$post2,
  ret = fates$ret,
  pil = fates$pil,
  study_week = (fates$study.week - mean(fates$study.week)) / sd(fates$study.week),
  
  # spline
  basis = basis
  
)

#_______________________________________________________________________________________________
# 4b. Parameters to monitor ----
#_______________________________________________________________________________________________

params <- c("a0",
            "a0_mean",
            "a0_sigma",
            "w_mean",
            "w_sigma",
            "w",
            "lambda",
            "hr_bci",
            "hr_bci_study_week",
            "hr_ret_total1", 
            "hr_ret_total2",
            "hr_pil_total1",
            "hr_pil_total2")

#_______________________________________________________________________________________________
# 4c. Initial values for MCMC ----
#_______________________________________________________________________________________________

# inits
inits <- list(
  
  # spline intercept
  a0 = rnorm(constants$nsf, log(0.02), 1), 
  a0_mean = rnorm(1, log(0.02), 1),   
  a0_sigma = rexp(1, 1),
  a0_z = rexp(constants$nsf, 1),
  
  # spline parameters (by sf:basis)
  w_mean = rnorm(1, 0, 1),
  w_sigma = rexp(1, 1),
  w_z = matrix(rnorm(constants$nsf * constants$nbasis, 0, 1),
               nrow = constants$nsf,
               ncol = constants$nbasis),
  w = matrix(rnorm(constants$nsf * constants$nbasis, 0, 1),
             nrow = constants$nsf,
             ncol = constants$nbasis),
  
  # smoothing penalty (by sf)
  lambda_raw_mean = rnorm(1, 0, 1),
  lambda_raw_sigma = rexp(1, 1),
  lambda_raw_z = rnorm(constants$nsf, 0, 1),
  
  # coefficients
  b_bci = rnorm(1, 0, 2.5),
  b_bci_study_week = rnorm(1, 0, 2.5),
  b_ret = rnorm(1, 0, 2.5),
  b_pil = rnorm(1, 0, 2.5),
  b_ret_post1 = rnorm(1, 0, 2.5),
  b_ret_post2 = rnorm(1, 0, 2.5),
  b_pil_post1 = rnorm(1, 0, 2.5),
  b_pil_post2 = rnorm(1, 0, 2.5)
  
)

#_______________________________________________________________________________________________
# 5. Model 1 (All morts, Scenario 1) ----
#_______________________________________________________________________________________________
# 5a. Code ----
#_______________________________________________________________________________________________

model.code.1 <- nimbleCode({
  
  # priors
  # baseline hazard spline
  # spline intercept
  a0_mean ~ dnorm(log(0.01), sd = 1)      # mean intercept
  a0_sigma ~ dexp(rate = 1)               # SD intercept
  
  # mean spline coefficient
  w_mean ~ dnorm(0, sd = 1)
  w_sigma ~ dexp(rate = 1)
  
  # smoothing penalty
  lambda_raw_mean ~ dnorm(0, sd = 1)
  lambda_raw_sigma ~ dexp(rate = 2)       # to help this sample, let's decrease this
  
  # coefficients (Cauchy priors)
  # intrinsic 
  b_bci ~ dt(0, sigma = 2.5, df = 1)
  
  # intrinsic / extrinsic interactions
  b_bci_study_week ~ dt(0, sigma = 2.5, df = 1)
  
  # treatment variables
  b_ret ~ dt(0, sigma = 2.5, df = 1)
  b_pil ~ dt(0, sigma = 2.5, df = 1)
  
  # pre-post treatment interaction variables
  b_ret_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_ret_post2 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post2 ~ dt(0, sigma = 2.5, df = 1)
  
  # non-centered random spline parameters (sex / forest type)
  for (y in 1:nsf) {
    
    # sex/forest type-specific intercepts
    a0_z[y] ~ dnorm(0, sd = 1)
    
    a0[y] <- a0_mean + a0_sigma * a0_z[y]
    
    # sex/forest type-specific smoothing penalty parameter (lambda)
    lambda_raw_z[y] ~ dnorm(0, sd = 1)
    
    lambda[y] <- exp(lambda_raw_mean + lambda_raw_sigma * lambda_raw_z[y])
    
    for (z in 1:nbasis) {
      
      w_z[y, z] ~ dnorm(0, sd = 1)
      
      w[y, z] <- w_mean + w_sigma * (w_z[y, z] / lambda[y])
      
    }
    
  }
  
  # likelihood
  for (i in 1:N) {
    
    # loop through basis functions because nimble is upset at non-scalar results
    for (j in 1:nbasis) {
      
      wb[i, j] <- w[s_f[i], j] * basis[i, j]
      
    }
    
    # all morts, scenario 1
    y_mort_1[i] ~ dpois(
      
      # baseline hazard (hierarchical spline)
      exp(
        
        a0[s_f[i]] + sum(wb[i, 1:nbasis])
        
        ) * 
        
      # covariate effects
      exp(
        
        b_bci * bci_1[i] +
        b_bci_study_week * bci_1[i] * study_week[i] +
        b_ret * ret[i] + 
        b_ret_post1 * post1[i] * ret[i] +
        b_ret_post2 * post2[i] * ret[i] +
        b_pil * pil[i] + 
        b_pil_post1 * post1[i] * pil[i] +
        b_pil_post2 * post2[i] * pil[i]
        
        )
      
      )
    
  }
  
  
  # derived quantities (hazard ratios)
  hr_bci <- exp(b_bci)
  hr_bci_study_week <- exp(b_bci_study_week)
  hr_ret_total1 <- exp(b_ret + b_ret_post1)
  hr_ret_total2 <- exp(b_ret + b_ret_post2)
  hr_pil_total1 <- exp(b_pil + b_pil_post1)
  hr_pil_total2 <- exp(b_pil + b_pil_post2)
  
})

#_______________________________________________________________________________________________
# 5b. Compile model ----
#_______________________________________________________________________________________________

# base model, with code, constants, data, and inits
model.1 <- nimbleModel(
  
  code = model.code.1,
  constants = constants,
  data = fates.1,
  inits = inits,
  buildDerivs = TRUE  # required for HMC
  
  ) 

# build HMC algorithm (include parameters to monitor)
model.hmc.1 <- buildHMC(model.1, monitors = params)

# compile initial model
model.comp.1 <- compileNimble(model.1)

# and now add in the HMC compilation
model.hmc.comp.1 <- compileNimble(model.hmc.1, project = model.1)

#_______________________________________________________________________________________________
# 5c. Run MCMC ----
#_______________________________________________________________________________________________

model.fit.1 <- runMCMC(
  
  mcmc = model.hmc.comp.1,          
  nchains = 3,                     
  nburnin = 7000,                   # 3000 total posterior samples
  niter = 8000, 
  samplesAsCodaMCMC = TRUE
  
  )

#_______________________________________________________________________________________________
# 5d. Save model output  ----
#_______________________________________________________________________________________________

saveRDS(model.fit.1, "Model outputs/model_1.rds")

#_______________________________________________________________________________________________
# 6. Model 2 (All morts, Scenario 2) ----
#_______________________________________________________________________________________________
# 6a. Code ----
#_______________________________________________________________________________________________

model.code.2 <- nimbleCode({
  
  # priors
  # baseline hazard spline
  # spline intercept
  a0_mean ~ dnorm(log(0.01), sd = 1)      # mean intercept
  a0_sigma ~ dexp(rate = 1)               # SD intercept
  
  # mean spline coefficient
  w_mean ~ dnorm(0, sd = 1)
  w_sigma ~ dexp(rate = 1)
  
  # smoothing penalty
  lambda_raw_mean ~ dnorm(0, sd = 1)
  lambda_raw_sigma ~ dexp(rate = 2)       # to help this sample, let's decrease this
  
  # coefficients (Cauchy priors)
  # intrinsic 
  b_bci ~ dt(0, sigma = 2.5, df = 1)
  
  # intrinsic / extrinsic interactions
  b_bci_study_week ~ dt(0, sigma = 2.5, df = 1)
  
  # treatment variables
  b_ret ~ dt(0, sigma = 2.5, df = 1)
  b_pil ~ dt(0, sigma = 2.5, df = 1)
  
  # pre-post treatment interaction variables
  b_ret_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_ret_post2 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post2 ~ dt(0, sigma = 2.5, df = 1)
  
  # non-centered random spline parameters (sex / forest type)
  for (y in 1:nsf) {
    
    # cluster-specific intercepts
    a0_z[y] ~ dnorm(0, sd = 1)
    
    a0[y] <- a0_mean + a0_sigma * a0_z[y]
    
    # cluster-specific smoothing penalty parameter (lambda)
    lambda_raw_z[y] ~ dnorm(0, sd = 1)
    
    lambda[y] <- exp(lambda_raw_mean + lambda_raw_sigma * lambda_raw_z[y])
    
    for (z in 1:nbasis) {
      
      w_z[y, z] ~ dnorm(0, sd = 1)
      
      w[y, z] <- w_mean + w_sigma * (w_z[y, z] / lambda[y])
      
    }
    
  }
  
  # likelihood
  for (i in 1:N) {
    
    # loop through basis functions because nimble is upset at non-scalar results
    for (j in 1:nbasis) {
      
      wb[i, j] <- w[s_f[i], j] * basis[i, j]
      
    }
    
    # all morts, scenario 2
    y_mort_2[i] ~ dpois(
      
      # baseline hazard (hierarchical spline)
      exp(
        
        a0[s_f[i]] + sum(wb[i, 1:nbasis])
        
      ) * 
        
        # covariate effects
        exp(
          
          b_bci * bci_1[i] +
            b_bci_study_week * bci_1[i] * study_week[i] +
            b_ret * ret[i] + 
            b_ret_post1 * post1[i] * ret[i] +
            b_ret_post2 * post2[i] * ret[i] +
            b_pil * pil[i] + 
            b_pil_post1 * post1[i] * pil[i] +
            b_pil_post2 * post2[i] * pil[i]
          
        )
      
    )
    
  }
  
  
  # derived quantities (hazard ratios)
  hr_bci <- exp(b_bci)
  hr_bci_study_week <- exp(b_bci_study_week)
  hr_ret_total1 <- exp(b_ret + b_ret_post1)
  hr_ret_total2 <- exp(b_ret + b_ret_post2)
  hr_pil_total1 <- exp(b_pil + b_pil_post1)
  hr_pil_total2 <- exp(b_pil + b_pil_post2)
  
})

#_______________________________________________________________________________________________
# 6b. Compile model ----
#_______________________________________________________________________________________________

# base model, with code, constants, data, and inits
model.2 <- nimbleModel(
  
  code = model.code.2,
  constants = constants,
  data = fates.1,
  inits = inits,
  buildDerivs = TRUE
  
  ) 

# build HMC algorithm (include parameters to monitor)
model.hmc.2 <- buildHMC(model.2, monitors = params)

# compile initial model
model.comp.2 <- compileNimble(model.2)

# and now add in the HMC compilation
model.hmc.comp.2 <- compileNimble(model.hmc.2, project = model.2)

#_______________________________________________________________________________________________
# 6c. Run MCMC ----
#_______________________________________________________________________________________________

model.fit.2 <- runMCMC(
  
  mcmc = model.hmc.comp.2,          
  nchains = 3,                     
  nburnin = 7000,                  
  niter = 8000,  
  samplesAsCodaMCMC = TRUE
  
)

#_______________________________________________________________________________________________
# 6d. Save model output  ----
#_______________________________________________________________________________________________

saveRDS(model.fit.2, "Model outputs/model_2.rds")

#_______________________________________________________________________________________________
# 7. Model 3 (All morts, Scenario 3) ----
#_______________________________________________________________________________________________
# 7a. Code ----
#_______________________________________________________________________________________________

model.code.3 <- nimbleCode({
  
  # priors
  # baseline hazard spline
  # spline intercept
  a0_mean ~ dnorm(log(0.01), sd = 1)      # mean intercept
  a0_sigma ~ dexp(rate = 1)               # SD intercept
  
  # mean spline coefficient
  w_mean ~ dnorm(0, sd = 1)
  w_sigma ~ dexp(rate = 1)
  
  # smoothing penalty
  lambda_raw_mean ~ dnorm(0, sd = 1)
  lambda_raw_sigma ~ dexp(rate = 2)       # to help this sample, let's decrease this
  
  # coefficients (Cauchy priors)
  # intrinsic 
  b_bci ~ dt(0, sigma = 2.5, df = 1)
  
  # intrinsic / extrinsic interactions
  b_bci_study_week ~ dt(0, sigma = 2.5, df = 1)
  
  # treatment variables
  b_ret ~ dt(0, sigma = 2.5, df = 1)
  b_pil ~ dt(0, sigma = 2.5, df = 1)
  
  # pre-post treatment interaction variables
  b_ret_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_ret_post2 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post2 ~ dt(0, sigma = 2.5, df = 1)
  
  # non-centered random spline parameters (sex / forest type)
  for (y in 1:nsf) {
    
    # cluster-specific intercepts
    a0_z[y] ~ dnorm(0, sd = 1)
    
    a0[y] <- a0_mean + a0_sigma * a0_z[y]
    
    # cluster-specific smoothing penalty parameter (lambda)
    lambda_raw_z[y] ~ dnorm(0, sd = 1)
    
    lambda[y] <- exp(lambda_raw_mean + lambda_raw_sigma * lambda_raw_z[y])
    
    for (z in 1:nbasis) {
      
      w_z[y, z] ~ dnorm(0, sd = 1)
      
      w[y, z] <- w_mean + w_sigma * (w_z[y, z] / lambda[y])
      
    }
    
  }
  
  # likelihood
  for (i in 1:N) {
    
    # loop through basis functions because nimble is upset at non-scalar results
    for (j in 1:nbasis) {
      
      wb[i, j] <- w[s_f[i], j] * basis[i, j]
      
    }
    
    # all morts, scenario 3
    y_mort_3[i] ~ dpois(
      
      # baseline hazard (hierarchical spline)
      exp(
        
        a0[s_f[i]] + sum(wb[i, 1:nbasis])
        
      ) * 
        
        # covariate effects
        exp(
          
          b_bci * bci_1[i] +
            b_bci_study_week * bci_1[i] * study_week[i] +
            b_ret * ret[i] + 
            b_ret_post1 * post1[i] * ret[i] +
            b_ret_post2 * post2[i] * ret[i] +
            b_pil * pil[i] + 
            b_pil_post1 * post1[i] * pil[i] +
            b_pil_post2 * post2[i] * pil[i]
          
        )
      
    )
    
  }
  
  
  # derived quantities (hazard ratios)
  hr_bci <- exp(b_bci)
  hr_bci_study_week <- exp(b_bci_study_week)
  hr_ret_total1 <- exp(b_ret + b_ret_post1)
  hr_ret_total2 <- exp(b_ret + b_ret_post2)
  hr_pil_total1 <- exp(b_pil + b_pil_post1)
  hr_pil_total2 <- exp(b_pil + b_pil_post2)
  
})

#_______________________________________________________________________________________________
# 7b. Compile model ----
#_______________________________________________________________________________________________

# base model, with code, constants, data, and inits
model.3 <- nimbleModel(
  
  code = model.code.3,
  constants = constants,
  data = fates.1,
  inits = inits,
  buildDerivs = TRUE
  
  )

# build HMC algorithm (include parameters to monitor)
model.hmc.3 <- buildHMC(model.3, monitors = params)

# compile initial model
model.comp.3 <- compileNimble(model.3)

# and now add in the HMC compilation
model.hmc.comp.3 <- compileNimble(model.hmc.3, project = model.3)

#_______________________________________________________________________________________________
# 7c. Run MCMC ----
#_______________________________________________________________________________________________

model.fit.3 <- runMCMC(
  
  mcmc = model.hmc.comp.3,          
  nchains = 3,                     
  nburnin = 7000,                 
  niter = 8000,          
  samplesAsCodaMCMC = TRUE
  
)

#_______________________________________________________________________________________________
# 7d. Save model output  ----
#_______________________________________________________________________________________________

saveRDS(model.fit.3, "Model outputs/model_3.rds")
