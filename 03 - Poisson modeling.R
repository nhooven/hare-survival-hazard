# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03b - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 29 Dec 2023
# Date last modified: 09 Jan 2024
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

fates <- read.csv("fates_cleaned_2.csv")

# keep only columns we need for modeling
fates.1 <- fates %>% dplyr::select(cluster,
                                   Site,
                                   Ear.tag,
                                   Collar.type,
                                   week,
                                   status.num,
                                   Sex.1,
                                   Mass.1,
                                   HFL.1,
                                   PC1,
                                   BCI.1,
                                   Treatment.Retention,
                                   Treatment.Piling)

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----
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
                     bsz = fates.1$PC1,
                     bci = fates.1$BCI,
                     ret = fates.1$Treatment.Retention,
                     pil = fates.1$Treatment.Piling,
                     clust = fates.1$cluster,
                     nclust = 4,
                     basis = basis,
                     num_basis = num_basis)

#_______________________________________________________________________________________________
# 4. M5 - Model with covariate effects and a spline on hazard ----
#_______________________________________________________________________________________________

m5 <- rstan::stan(
  
  model_code = "
  
  data {
  
  // number of observations
  int N;                       // obs
  int num_basis;               // basis functions
  int nclust;                  // clusters
  
  // response variable (binary, 0-1)
  int y[N];
  
  // time (required for estimating time-variant hazard);
  // time in week of the year (1-52)
  vector[N] t;       
  
  // spline matrix
  matrix[N, num_basis] basis;
  
  // covariates (continuous)
  real bsz[N];             // body size principal component (unitless) 
  real bci[N];             // standardized body condition index (mass/hfl)
  
  // covariates (categorical)
  int clust[N];            // index for cluster (1-4)
  int sex[N];              // indicator for sex (0 = F, 1 = M)
  int ret[N];              // indicator for retention
  int pil[N];              // indicator for piling
  
  }
  
  parameters {
  
  real a0;                          // spline intercept
  vector[num_basis] w0;             // weights
  vector[nclust] c0;                // intercept by cluster
  real b_sex;                       // slope for sex
  real b_ret;                       // slope for retention
  real b_pil;                       // slope for piling
  real b_bsz;                       // slope for body size
  real b_bci;                       // slope for body condition
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  w0 ~ normal(0, 1);                     // weights
  a0 ~ normal(0, 1);                     // intercept
  
  // coefficients
  c0[clust] ~ normal(0, 1);               // normal prior on c0
  b_sex ~ normal(0, 2.5);                // normal prior on b_sex
  b_ret ~ normal(0, 2.5);                // normal prior on b_ret
  b_pil ~ normal(0, 2.5);                // normal prior on b_pil
  b_bsz ~ normal(0, 2.5);                // normal prior on b_bsz
  b_bci ~ normal(0, 2.5);                // normal prior on b_bci
  
  // model
  // linear predictor
  vector[N] bw;
  vector[N] lambda;
  
  bw = to_vector(basis * w0);
  
  for (i in 1:N) {
  
  lambda[i] = exp(a0 + bw[i]) *
  
            exp(c0[clust[i]] +
                b_sex * sex[i] + 
                b_ret * ret[i] + 
                b_pil * pil[i] +
                b_bsz * bsz[i] +
                b_bci * bci[i]);
  
  }
  
  y ~ poisson(lambda);
  
  }
  
  generated quantities {
  
  // hazard ratios
  real hr_sex = exp(b_sex);
  real hr_ret = exp(b_ret);
  real hr_pil = exp(b_pil);
  real hr_bsz = exp(b_bsz);
  real hr_bci = exp(b_bci);
  
  }
  
  ",
  
  data = fates.stan.1,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m5)

# estimates
plot(m5, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_bsz", "hr_bci"))
plot(m5, pars = c("c0[1]", "c0[2]", "c0[3]", "c0[4]"))

# trace
traceplot(m5, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_bsz", "hr_bci"))
traceplot(m5, pars = c("c0[1]", "c0[2]", "c0[3]", "c0[4]"))

#_______________________________________________________________________________________________
# 5. Save image ----
#_______________________________________________________________________________________________

save.image("poisson_model.RData")
