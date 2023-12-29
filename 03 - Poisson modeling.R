# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 29 Dec 2023
# Date last modified: 29 Dec 2023
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
  y ~ poisson(exp(a0 + to_vector(basis * w0)));                           
  
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

#_______________________________________________________________________________________________
# 7. Examine the spline predictions ----
#_______________________________________________________________________________________________

# extract draws
m4.draws <- rstan::extract(m4)

# first, we need to create the "normal" scale predictions
# we'll start by using the posterior medians and 95% quantiles, then loop through all of them (apply?)
m4.50 <- data.frame(a0 = median(m4.draws$a0),
                    w0.1 = median(m4.draws$w0[ ,1]),
                    w0.2 = median(m4.draws$w0[ ,2]),
                    w0.3 = median(m4.draws$w0[ ,3]),
                    w0.4 = median(m4.draws$w0[ ,4]),
                    w0.5 = median(m4.draws$w0[ ,5]),
                    w0.6 = median(m4.draws$w0[ ,6]))

m4.lo <- data.frame(a0 = quantile(m4.draws$a0, prob = 0.025),
                    w0.1 = quantile(m4.draws$w0[ ,1], prob = 0.25),
                    w0.2 = quantile(m4.draws$w0[ ,2], prob = 0.25),
                    w0.3 = quantile(m4.draws$w0[ ,3], prob = 0.25),
                    w0.4 = quantile(m4.draws$w0[ ,4], prob = 0.25),
                    w0.5 = quantile(m4.draws$w0[ ,5], prob = 0.25),
                    w0.6 = quantile(m4.draws$w0[ ,6], prob = 0.25))

m4.up <- data.frame(a0 = quantile(m4.draws$a0, prob = 0.975),
                    w0.1 = quantile(m4.draws$w0[ ,1], prob = 0.75),
                    w0.2 = quantile(m4.draws$w0[ ,2], prob = 0.75),
                    w0.3 = quantile(m4.draws$w0[ ,3], prob = 0.75),
                    w0.4 = quantile(m4.draws$w0[ ,4], prob = 0.75),
                    w0.5 = quantile(m4.draws$w0[ ,5], prob = 0.75),
                    w0.6 = quantile(m4.draws$w0[ ,6], prob = 0.75))

# now multiply the "normal" scale weight vector by the basis function matrix
m4.w0.vec.50 <- t(as.matrix(as.numeric(m4.50[ ,2:7])))
m4.w0.vec.lo <- t(as.matrix(as.numeric(m4.lo[ ,2:7])))
m4.w0.vec.up <- t(as.matrix(as.numeric(m4.up[ ,2:7])))

w.by.basis.50 <- sweep(basis, 2, m4.w0.vec.50, `*`)
w.by.basis.lo <- sweep(basis, 2, m4.w0.vec.lo, `*`)
w.by.basis.up <- sweep(basis, 2, m4.w0.vec.up, `*`)

w.by.basis.sum.50 <- apply(w.by.basis.50, 1, sum)
w.by.basis.sum.lo <- apply(w.by.basis.lo, 1, sum)
w.by.basis.sum.up <- apply(w.by.basis.up, 1, sum)

# add to the intercept (a0)
spline.pred.50 <- m4.50$a0 + w.by.basis.sum.50
spline.pred.lo <- m4.lo$a0 + w.by.basis.sum.lo
spline.pred.up <- m4.up$a0 + w.by.basis.sum.up

# add to data.frame and exponentiate (to transform to hazard scale)
pred.data <- data.frame(y = exp(spline.pred.50),
                        lo = exp(spline.pred.lo),
                        up = exp(spline.pred.up),
                        week = fates.2$week)

# change fates to max of the spline
fates.pred <- fates.2 %>% mutate(mort = ifelse(status.num == 1,
                                               max(pred.data$up),
                                               min(pred.data$lo)))

ggplot() +
  
  # white background
  theme_bw() +
  
  # line for spline prediction
  geom_line(data = pred.data,
            aes(x = week,
                y = y),     
            linewidth = 1.25) +
  
  # ribbon for credible intervals
  geom_ribbon(data = pred.data,
            aes(x = week,
                y = y,
                ymin = lo,
                ymax = up),     
            linewidth = 1.25,
            alpha = 0.25) +
  
  # add jittered points for morts/non-morts
  geom_jitter(data = fates.pred,
             aes(x = week,
                 y = mort),
             height = 0.0001,
             alpha = 0.05) +
  
  # axis labels
  ylab("Baseline hazard") +
  xlab("Week of year")
  
#_______________________________________________________________________________________________
# 8. M5 - Model with covariate effects and a spline on hazard ----
#_______________________________________________________________________________________________

m5 <- rstan::stan(
  
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
  
  // covariates (continuous)
  real mass[N];            // standardized mass (kg)
  real hfl[N];             // standardized hind foot length (cm)
  
  // covariates (categorical)
  int sex[N];              // indicator for sex (0 = F, 1 = M)
  int ret[N];              // indicator for retention
  int pil[N];              // indicator for piling
  
  }
  
  parameters {
  
  real a0;                          // intercept
  vector[num_basis] w0;             // weights
  real b_sex;                       // slope for sex
  real b_ret;                       // slope for retention
  real b_pil;                       // slope for piling
  real b_mas;                       // slope for mass
  real b_hfl;                       // slope for hfl
  
  }
  
  model {
  
  // priors (these are on a normal scale)
  // normal scale priors for spline parameters
  w0 ~ normal(0, 1);                     // weights
  a0 ~ normal(0, 1);                     // intercept
  
  // coefficients
  b_sex ~ normal(0, 2.5);                // normal prior on b_sex
  b_ret ~ normal(0, 2.5);                // normal prior on b_ret
  b_pil ~ normal(0, 2.5);                // normal prior on b_pil
  b_mas ~ normal(0, 2.5);                // normal prior on b_mas
  b_hfl ~ normal(0, 2.5);                // normal prior on b_hfl
  
  // model
  // linear predictor
  vector[N] bw;
  vector[N] lambda;
  
  bw = to_vector(basis * w0);
  
  for (i in 1:N) {
  
  lambda[i] = exp(a0 + bw[i]) *
  
            exp(b_sex * sex[i] + 
                b_ret * ret[i] + 
                b_pil * pil[i] +
                b_mas * mass[i] +
                b_hfl * hfl[i]);
  
  }
  
  y ~ poisson(lambda);
  
  }
  
  generated quantities {
  
  // hazard ratios
  real hr_sex = exp(b_sex);
  real hr_ret = exp(b_ret);
  real hr_pil = exp(b_pil);
  real hr_mas = exp(b_mas);
  real hr_hfl = exp(b_hfl);
  
  }
  
  ",
  
  data = fates.stan.1,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m5)
plot(m5, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_mas", "hr_hfl"))
traceplot(m5, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_mas", "hr_hfl"))

#_______________________________________________________________________________________________
# 9. Save image ----
#_______________________________________________________________________________________________

save.image("poisson_model.RData")
