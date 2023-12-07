# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("fates_cleaned.csv")

# keep only columns we need for modeling
fates.1 <- fates %>% dplyr::select(Site,
                                   Ear.tag,
                                   Collar.type,
                                   week,
                                   status.num,
                                   Sex.1,
                                   Mass.1,
                                   HFL.1,
                                   Treatment.1)

# data list for Stan
fates.stan <- list(N = nrow(fates.1),
                   y = fates.1$status.num,
                   t = fates.1$week,
                   sex = fates.1$Sex.1,
                   sex_lv = unique(fates.1$Sex.1),
                   mass = fates.1$Mass.1,
                   hfl = fates.1$HFL.1,
                   trt = fates.1$Treatment.1)

#_______________________________________________________________________________________________
# 3. Build models using Stan ----
#_______________________________________________________________________________________________

# Stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#_______________________________________________________________________________________________
# 3a. M1 - constant baseline hazard, no covariate effects (SANITY CHECK) ----
#_______________________________________________________________________________________________

m1 <- stan(
  
  model_code = "
  
  data {
  
  // number of observations
  int N;            
  
  // response variable (binary, 0-1)
  int y[N];                
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  int t[N];       
  
  // covariates
  int sex[N];              // index for sex (1 = F, 0 = M)
  real mass[N];            // standardized mass (kg)
  real hfl[N];             // standardized hind foot length (cm)
  int trt[N];              // index for treatment (1 = control, 2 = retention, 3 = piling)
  
  }
  
  parameters {
  
  real lambda0;            // baseline hazard
  
  }
  
  transformed parameters {
  
  real lambda;             // mean of the Poisson likelihood
  lambda = exp(lambda0);   // exponentiated baseline hazard
  
  }
  
  model {
  
  // priors
  target += gamma_lpdf(lambda | 1, 1);
  
  // Poisson likelihood
  y ~ poisson(lambda);
  
  
  }
  
  
  ",
  
  data = fates.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m1)
plot(m1)
traceplot(m1)

#_______________________________________________________________________________________________
# 3b. M2 - constant baseline hazard, sex effect ----
#_______________________________________________________________________________________________

# this was the first model that worked beautifully! I'll keep it here

m2 <- rstan::stan(
  
  model_code = "
  
  data {
  
  // number of observations
  int N;            
  
  // response variable (binary, 0-1)
  int y[N];                
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  int t[N];       
  
  // covariates
  real mass[N];            // standardized mass (kg)
  real hfl[N];             // standardized hind foot length (cm)
  int trt[N];              // index for treatment (1 = control, 2 = retention, 3 = piling)
  
  // sex (categorical)
  array[N] int sex;              // index for sex (1 = F, 2 = M)
  
  }
  
  parameters {
  
  real<lower = 0> h0;               // baseline hazard
  vector[2] b_sex;                  // slope for sex
  
  }
  
  model {
  
  // lambda
  vector[N] lambda;
  
  // priors
  h0 ~ gamma(2, 20);                     // gamma prior on baseline hazard
  b_sex ~ normal(0, 2.5);                // normal prior on b.sex

  // Poisson likelihood
  // Linear predictor
  for (i in 1:N) {
  
        lambda[i] = exp(h0 + b_sex[sex[i]] * sex[i]);
        
  }
    
  y ~ poisson(lambda);
  
  }
  
  ",
  
  data = fates.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m2)
plot(m2)
traceplot(m2)

#_______________________________________________________________________________________________
# 3c. M3 - constant baseline hazard, all covariates ----
#_______________________________________________________________________________________________

# this was the first model that worked beautifully! I'll keep it here

m3 <- rstan::stan(
  
  model_code = "
  
  data {
  
  // number of observations
  int N;            
  
  // response variable (binary, 0-1)
  int y[N];                
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  int t[N];       
  
  // covariates (continuous)
  real mass[N];            // standardized mass (kg)
  real hfl[N];             // standardized hind foot length (cm)
  
  // covariates (categorical)
  array[N] int sex;              // index for sex (1 = F, 2 = M)
  array[N] int trt;              // index for treatment (1 = control, 2 = retention, 3 = piling)
  
  }
  
  parameters {
  
  real<lower = 0> h0;               // baseline hazard
  vector[2] b_sex;                  // slope for sex
  vector[3] b_trt;                  // slope for treatment
  real b_mas;                       // slope for mass
  real b_hfl;                       // slope for hfl
  
  }
  
  model {
  
  // lambda
  vector[N] lambda;
  
  // priors
  h0 ~ gamma(2, 20);                     // gamma prior on baseline hazard
  b_sex ~ normal(0, 2.5);                // normal prior on b_sex
  b_trt ~ normal(0, 2.5);                // normal prior on b_trt
  b_mas ~ normal(0, 2.5);                // normal prior on b_mas
  b_hfl ~ normal(0, 2.5);                // normal prior on b_hfl

  // Poisson likelihood
  // Linear predictor
  for (i in 1:N) {
  
        lambda[i] = exp(h0 + 
                        b_sex[sex[i]] * sex[i] + 
                        b_trt[trt[i]] * trt[i] + 
                        b_mas * mass[i] +
                        b_hfl * hfl[i]);
        
  }
    
  y ~ poisson(lambda);
  
  }
  
  ",
  
  data = fates.stan,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m3)
plot(m3)
traceplot(m3)

# these look pretty good - interesting that the sex effects dampen with the inclusion of the
# other predictors

# I suspect that there is some multicollinearity going on here
# also - should I change index to indicator variables here?
