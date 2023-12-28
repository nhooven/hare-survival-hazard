# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 
# Date last modified: 28 Dec 2023
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)         # construct basic functions

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
                                   Treatment.Retention,
                                   Treatment.Piling)

# data list for Stan
fates.stan <- list(N = nrow(fates.1),
                   y = fates.1$status.num,
                   t = fates.1$week,
                   sex = fates.1$Sex.1,
                   mass = fates.1$Mass.1,
                   hfl = fates.1$HFL.1,
                   ret = fates.1$Treatment.Retention,
                   pil = fates.1$Treatment.Piling)

#_______________________________________________________________________________________________
# 3. Examine potential multicollinearity ----
#_______________________________________________________________________________________________

# correlation between continuous covariates (mass and HFL)
cor(x = fates.1[ ,c("Mass.1", "HFL.1")], method = "pearson")

ggplot(data = fates.1,
       aes(x = Mass.1,
           y = HFL.1)) +
  theme_bw() +
  geom_point()

# kruskal-wallis tests for continuous by categorical
# mass by sex
kruskal.test(x = fates.1$Mass.1, g = fates.1$Sex.1)

ggplot(data = fates.1,
       aes(fill = as.factor(Sex.1),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by sex
kruskal.test(x = fates.1$HFL.1, g = fates.1$Sex.1)

ggplot(data = fates.1,
       aes(fill = as.factor(Sex.1),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# mass by retention
kruskal.test(x = fates.1$Mass.1, g = fates.1$Treatment.Retention)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Retention),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by retention
kruskal.test(x = fates.1$HFL.1, g = fates.1$Treatment.Retention)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Retention),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# mass by piling
kruskal.test(x = fates.1$Mass.1, g = fates.1$Treatment.Piling)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Piling),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by Piling
kruskal.test(x = fates.1$HFL.1, g = fates.1$Treatment.Piling)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Piling),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

#_______________________________________________________________________________________________
# 4. Build models using Stan ----
#_______________________________________________________________________________________________

# Stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#_______________________________________________________________________________________________
# 4a. M1 - constant baseline hazard, no covariate effects (SANITY CHECK) ----
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
# 4b. M2 - constant baseline hazard, sex effect ----
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
# 4c. M3 - constant baseline hazard, all covariates ----
#_______________________________________________________________________________________________

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
  array[N] int sex;              // indicator for sex (0 = F, 1 = M)
  array[N] int ret;              // indicator for retention
  array[N] int pil;              // indicator for piling
  
  }
  
  parameters {
  
  real<lower = 0> h0;               // baseline hazard
  real b_sex;                       // slope for sex
  real b_ret;                       // slope for retention
  real b_pil;                       // slope for piling
  real b_mas;                       // slope for mass
  real b_hfl;                       // slope for hfl
  
  }
  
  model {
  
  // lambda
  vector[N] lambda;
  
  // priors
  h0 ~ gamma(2, 20);                     // gamma prior on baseline hazard
  b_sex ~ normal(0, 2.5);                // normal prior on b_sex
  b_ret ~ normal(0, 2.5);                // normal prior on b_ret
  b_pil ~ normal(0, 2.5);                // normal prior on b_pil
  b_mas ~ normal(0, 2.5);                // normal prior on b_mas
  b_hfl ~ normal(0, 2.5);                // normal prior on b_hfl

  // Poisson likelihood
  // Linear predictor
  for (i in 1:N) {
  
        lambda[i] = h0 * exp(b_sex * sex[i] + 
                             b_ret * ret[i] + 
                             b_pil * pil[i] +
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

# this looks good - let's include a spline on the baseline hazard!

#_______________________________________________________________________________________________
# 5. Construct basis functions for spline on hazard ----
#_______________________________________________________________________________________________

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

# define knots (quantile)
n.knots <- 5

knot.list <- quantile(fates.1$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# creating the basis matrix
basis <- bs(fates.1$week,                   # create weeks 
            knots = knot.list[-c(1, n.knots)],          
            degree = 3,                     # cubic spline
            intercept = FALSE)

# plot basis functions
basis.plot <- tibble(week = fates.1$week,
                     b1 = basis[ , 1],
                     b2 = basis[ , 2],
                     b3 = basis[ , 3],
                     b4 = basis[ , 4],
                     b5 = basis[ , 5],
                     b6 = basis[ , 6]) %>%
              pivot_longer(cols = b1:b6)

ggplot(data = fates.1,
       aes(x = week,
           y = status.num)) +
  theme_bw() +
  geom_jitter(alpha = 0.1,
              height = 0.01) +
  geom_line(data = basis.plot,
            aes(x = week,
                y = value,
                group = name),
            linewidth = 1.5,
            alpha = 0.25)

num_basis <- ncol(basis)

# add to new list
fates.stan.1 <- list(N = nrow(fates.1),
                     y = fates.1$status.num,
                     t = fates.1$week,
                     sex = fates.1$Sex.1,
                     mass = fates.1$Mass.1,
                     hfl = fates.1$HFL.1,
                     ret = fates.1$Treatment.Retention,
                     pil = fates.1$Treatment.Piling,
                     basis = basis,
                     num_basis = num_basis)

#_______________________________________________________________________________________________
# 6. M4 - Model with only the basis spline on hazard ----
#_______________________________________________________________________________________________

m4 <- rstan::stan(
  
  model_code = "
  
  data {
  
  // number of observations
  int N;
  int num_basis;
  
  // response variable (binary, 0-1)
  int y[N];
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  vector[N] t;       
  
  // spline matrix
  matrix[N, num_basis] basis;
  
  }
  
  parameters {
  
  real a0;                          // intercept
  vector[num_basis] w0;              // weights
  
  }
  
  transformed parameters {
  
  real h0;
  h0 = exp(a0);
  
  vector[num_basis] w;
  w = exp(w0);
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  w0 ~ normal(0, 1);                     // weights
  a0 ~ normal(0, 1);                     // intercept                  
  
  // model
  y ~ poisson(exp(a0 + to_vector(basis * w)));                           
  
  }",
  
  data = fates.stan.1,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m4)
plot(m4)
traceplot(m4)


# looks like i solved the sampling issues by:
# 1: constraining the hazard rates to be positive
# 2: transforming them to be positive by exponentiating

# now I need to calculate the entire spline first,
# then exponentiate it to find the spline on the scale we need
