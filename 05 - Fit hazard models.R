# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 05 - Fit hazard models
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 25 Nov 2025
# Date last modified: 25 Feb 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(nimble)          # modeling with nimble
library(nimbleHMC)       # HMC in nimble
library(coda)            # model outputs

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

constant.list <- readRDS("Cleaned data/constants.rds")
data.list <- readRDS("Cleaned data/data.rds")

# clone data.lists with different y to avoid having to write three identical models
data.list.1 <- data.list
data.list.2 <- data.list
data.list.3 <- data.list

data.list.1$y <- data.list$y_mort_1
data.list.2$y <- data.list$y_mort_2
data.list.3$y <- data.list$y_mort_3

#_______________________________________________________________________________________________
# 3. Model code ----
#_______________________________________________________________________________________________

model.code <- nimbleCode({
  
  # priors
  # baseline hazard spline
  # a0 - spline intercept
  a0_mean ~ dnorm(log(0.01), sd = 1)      # mean intercept
  a0_sigma ~ dexp(rate = 1)               # SD intercept
  
  # mean spline coefficient
  w_mean ~ dnorm(0, sd = 1)
  w_sigma ~ dexp(rate = 1)
  
  # lambda - smoothing penalty (log scale)
  loglam_mean ~ dnorm(0, sd = 1)
  loglam_sigma ~ dexp(rate = 2)       # to help this sample, let's decrease this
  
  # coefficients (Cauchy priors)
  # intrinsic 
  b_bci ~ dt(0, sigma = 2.5, df = 1)
  
  # intrinsic / extrinsic interactions
  b_bci_study_week ~ dt(0, sigma = 2.5, df = 1)
  
  # landscape
  b_dm ~ dt(0, sigma = 2.5, df = 1)
  b_open ~ dt(0, sigma = 2.5, df = 1)
  
  # treatment variables
  b_ret ~ dt(0, sigma = 2.5, df = 1)
  b_pil ~ dt(0, sigma = 2.5, df = 1)
  
  # pre-post treatment interaction variables
  b_ret_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_ret_post2 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post1 ~ dt(0, sigma = 2.5, df = 1)
  b_pil_post2 ~ dt(0, sigma = 2.5, df = 1)
  
  # random spline parameters - cluster (c) nested within sex (s)
  for (s in 1:2) {
    
    # parameters by sex
    # a0s - intercepts [2]
    a0s_z[s] ~ dnorm(0, sd = 1)
    
    a0s[s] <- a0_mean + a0_sigma * a0s_z[s]
    
    # lams - lambda - smoothing penalties [2]
    lams_z[s] ~ dnorm(0, sd = 1)
    
    log(lams[s]) <- loglam_mean + loglam_sigma * lams_z[s]
    
    # ws - spline weights [2 x nbasis]
    for (b in 1:nbasis) {
      
      ws_z[s, b] ~ dnorm(0, sd = 1)
      
      ws[s, b] <- w_mean + w_sigma * (ws_z[s, b] / lams[s])
      
    } # basis
    
    for (c in 1:nclust) {
      
      # hyperparameters by cluster, nested within sex
      # a0sc - intercepts [2 x 4]
      a0sc[s, c] ~ dnorm(a0s[s], sd = 1)
      
      # lamsc - lambda - smoothing penalties [2 x 4]
      lamsc[s, c] ~ dnorm(log(lams[s]), sd = 1)
      
      # wsc - spline weights [2 x 4 x nbasis]
      for (b in 1:nbasis) {
        
        wsc_z[s, c, b] ~ dnorm(ws_z[s, b], sd = 1)
        
        wsc[s, c, b] <- ws[s, b] + (wsc_z[s, c, b] / lamsc[s, c]) # assume sigma is 1 here
        
      } # basis
      
    } # cluster
    
  } # sex
  
  # likelihood calculations
  for (i in 1:N) {
    
    # weights x basis functions
    wb[i, 1:nbasis] <- wsc[sex[i], cluster[i], 1:nbasis] * basis[i, 1:nbasis]
    
    # Poisson likelihood
    y[i] ~ dpois(
      
      # baseline hazard (hierarchical spline)
      exp(a0sc[sex[i], cluster[i]] + sum(wb[i, 1:nbasis])) * 
        
        # covariate effects
        exp(
          
          b_bci * bci_1[i] +
          b_bci_study_week * bci_1[i] * study_week[i] +
          b_ret * ret[i] + 
          b_ret_post1 * post1[i] * ret[i] +
          b_ret_post2 * post2[i] * ret[i] +
          b_pil * pil[i] + 
          b_pil_post1 * post1[i] * pil[i] +
          b_pil_post2 * post2[i] * pil[i] +
          b_dm * dm[i] +
          b_open * open[i]
          
        )
      
    )
    
  } # N
  
  # derived quantities (hazard ratios)
  hr_bci <- exp(b_bci)
  hr_bci_study_week <- exp(b_bci_study_week)
  hr_ret_total1 <- exp(b_ret + b_ret_post1)
  hr_ret_total2 <- exp(b_ret + b_ret_post2)
  hr_pil_total1 <- exp(b_pil + b_pil_post1)
  hr_pil_total2 <- exp(b_pil + b_pil_post2)
  hr_dm <- exp(b_dm)
  hr_open <- exp(b_open)
    
    
})

# ______________________________________________________________________________
# 4. Initial values ----
# ______________________________________________________________________________

inits <- list(
  
  # baseline hazard spline
  # a0 - spline intercept
  a0_mean = rnorm(1, log(0.01), 1),      # mean intercept
  a0_sigma = rexp(1, 1),                # SD intercept
  
  # w - mean spline coefficient
  w_mean = rnorm(1, 0, sd = 1),
  w_sigma = rexp(1, rate = 1),
  
  # lambda - smoothing penalty (log scale)
  loglam_mean = rnorm(1, 0, sd = 1),
  loglam_sigma = rexp(1, rate = 2),       
  
  # coefficients (Cauchy priors)
  b_bci = rnorm(1, 0, 2.5),
  b_bci_study_week = rnorm(1, 0, 2.5),
  b_dm = rnorm(1, 0, 2.5),
  b_open = rnorm(1, 0, 2.5),
  b_ret = rnorm(1, 0, 2.5),
  b_pil = rnorm(1, 0, 2.5),
  b_ret_post1 = rnorm(1, 0, 2.5),
  b_ret_post2 = rnorm(1, 0, 2.5),
  b_pil_post1 = rnorm(1, 0, 2.5),
  b_pil_post2 = rnorm(1, 0, 2.5),
  
  # random spline parameters
  # sex
  a0s_z = rnorm(2, 0, sd = 1),
  lams_z = rnorm(2, 0, sd = 1),
    
  # sex x basis [2, nbasis]
  ws_z = matrix(rnorm(2 * constant.list$nbasis, 0, 1),
                nrow = 2,
                ncol = constant.list$nbasis),

  # sex / cluster [2, 4]
  a0sc = matrix(rnorm(2 * 4, 0, sd = 1),
                nrow = 2,
                ncol = 4),
  lamsc = matrix(rnorm(2 * 4, 0, sd = 1),
                 nrow = 2,
                 ncol = 4),
  
  # (sex / cluster) x basis [2, 4, nbasis]
  wsc_z = array(rnorm(2 * 4 * constant.list$nbasis, 0, sd = 1),
                dim = c(2, 4, constant.list$nbasis))
  
)

# ______________________________________________________________________________
# 5. Parameters to monitor ----
# ______________________________________________________________________________

monitor <- c(
  
  # mean spline parameters
  "a0_mean",
  "a0_sigma",
  "w_mean",
  "w_sigma",
  
  # random spline parameters
  "a0s",
  "ws",
  "a0sc",
  "wsc",
  
  # hazard ratios
  "hr_bci",
  "hr_bci_study_week",
  "hr_ret_total1", 
  "hr_ret_total2",
  "hr_pil_total1",
  "hr_pil_total2",
  "hr_dm",
  "hr_open"
  
  )

# ______________________________________________________________________________
# 6. Set up models ----
# ______________________________________________________________________________

# model 1
model.1 <- nimbleModel(
  
  code = model.code,
  constants = constant.list,
  data = data.list.1,
  inits = inits,
  buildDerivs = TRUE  
  
)

# model 2
model.2 <- nimbleModel(
  
  code = model.code,
  constants = constant.list,
  data = data.list.2,
  inits = inits,
  buildDerivs = TRUE  
  
)

# model 3
model.3 <- nimbleModel(
  
  code = model.code,
  constants = constant.list,
  data = data.list.3,
  inits = inits,
  buildDerivs = TRUE  
  
)

# check node size
length(model.1$getNodeNames())

# ______________________________________________________________________________
# 7. Compile models ----
# ______________________________________________________________________________

model.1.comp <- compileNimble(model.1)
model.2.comp <- compileNimble(model.2)
model.3.comp <- compileNimble(model.3)

# ______________________________________________________________________________
# 8. Build MCMC ----

# must include parameters to monitor

# ______________________________________________________________________________

model.1.mcmc <- buildHMC(model.1, monitors = monitor)
model.2.mcmc <- buildHMC(model.2, monitors = monitor)
model.3.mcmc <- buildHMC(model.3, monitors = monitor)

# ______________________________________________________________________________
# 9. Compile MCMC ----
# ______________________________________________________________________________

mcmc.1.comp <- compileNimble(model.1.mcmc, project = model.1)
mcmc.2.comp <- compileNimble(model.2.mcmc, project = model.2)
mcmc.3.comp <- compileNimble(model.3.mcmc, project = model.3)

# ______________________________________________________________________________
# 10. Run sampling ----

# we'll aim for 3,000 posterior samples, keeping 1k from each chain

# ______________________________________________________________________________

# model 1
model.1.run <- runMCMC(
  
  mcmc = mcmc.1.comp,          
  nchains = 3,                     
  nburnin = 7000,                   
  niter = 8000, 
  samplesAsCodaMCMC = TRUE
  
)

# model 2
model.2.run <- runMCMC(
  
  mcmc = mcmc.2.comp,          
  nchains = 3,                     
  nburnin = 7000,                   
  niter = 8000, 
  samplesAsCodaMCMC = TRUE
  
)

# model 3
model.3.run <- runMCMC(
  
  mcmc = mcmc.3.comp,          
  nchains = 3,                     
  nburnin = 7000,                   
  niter = 8000, 
  samplesAsCodaMCMC = TRUE
  
)

# ______________________________________________________________________________
# 11. Save samples ----
# ______________________________________________________________________________

saveRDS(model.1.run, "Model outputs/model_1.rds")
saveRDS(model.2.run, "Model outputs/model_2.rds")
saveRDS(model.3.run, "Model outputs/model_3.rds")
