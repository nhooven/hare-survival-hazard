# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 05 - Posterior predictive checks
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 05 Jan 2024
# Date completed: 
# Date last modified: 05 Jan 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

load("poisson_model.RData")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(m5))

#_______________________________________________________________________________________________
# 4. Distributions of lifetimes (in weeks) ----
#_______________________________________________________________________________________________
# 4a. Empirical distribution ----
#_______________________________________________________________________________________________

# loop through each individual, calculate lifetime (in weeks), associate covariates as needed
lifetimes.empirical <- data.frame()

for (i in 1:length(unique(fates.1$Ear.tag))) {
  
  # subset
  focal.indiv <- unique(fates.1$Ear.tag)[i]
  focal.df <- fates.1 %>% filter(Ear.tag == focal.indiv)
  
  # df of individual information
  lifetime.indiv <- data.frame(indiv = focal.indiv,
                               enter = focal.df$week[1],
                               exit = focal.df$week[nrow(focal.df)],
                               lifetime = nrow(focal.df),
                               died = ifelse(sum(focal.df$status.num) == 1,
                                             1,
                                             0),
                               sex = focal.df$Sex.1[1],
                               mass = focal.df$Mass.1[1])
  
  # bind into master df
  lifetimes.empirical <- rbind(lifetimes.empirical, lifetime.indiv)
  
}

# all monitoring times
ggplot(lifetimes.empirical,
       aes(x = lifetime)) +
  
  theme_bw() +
  
  geom_density(linewidth = 1.25)

# monitoring times - split by fate
ggplot(lifetimes.empirical,
       aes(x = lifetime,
           group = as.factor(died),
           fill = as.factor(died))) +
  
  theme_bw() +
  
  geom_density(linewidth = 1,
               alpha = 0.5)

# lifetimes related to mass
ggplot(lifetimes.empirical %>% filter(died == 1),
       aes(x = mass,
           y = lifetime)) +
  
  theme_bw() +
  
  geom_point()

# enter times (we'll draw the enter times from the simulation from this distribution)
ggplot(lifetimes.empirical,
       aes(x = enter)) +
  
  theme_bw() +
  
  geom_histogram(binwidth = 1,
                 color = "black",
                 fill = NA)

# now, let's look at the censoring rate (probability that the last observation is a zero)
# crudely, this will be the number of censor events * total weeks
censor.rate <- (nrow(lifetimes.empirical) - 
                sum(lifetimes.empirical$died)) / 
  sum(lifetimes.empirical$lifetime)

censor.rate

# thus the "hazard rate" of being censored is ~ 0.04
censor.prob <- exp(-censor.rate)

censor.prob

#_______________________________________________________________________________________________
# 4b. Simulate distributions of individuals from the posterior ----

# define function that starts counting from 1 again at the final step (typically week 52)
recycle.time <- function(enter = 1,
                         max = n.steps) {
  
  # define vector of length n.steps to hold 
  vec <- vector(mode = "integer", length = n.steps)
  
  # add enter time
  vec[1] <- enter
  
  # add 1
  for (i in 1:length(vec)) {
    
    # don't change the 1st entry
    if (i > 1) {
      
      # add one to all places to the right of the entry
      vec[i] <- vec[i - 1] + 1
      
      # subtract the n.steps from anything greater than n.steps
      vec[i] <- ifelse(vec[i] > n.steps,
                       vec[i] - n.steps,
                       vec[i])
      
    }
    
  }
  
  return(vec)
  
}

#_______________________________________________________________________________________________

# number of simulations, individuals, and simulation timesteps
n.sim <- 10
n.obs <- 200
n.steps <- 52

# let's create a prediction df first
pred.df <- data.frame(id = 1:n.obs)

# now let's draw enter times and covariate values from empirical distributions
set.seed(1234)

pred.df <- pred.df %>% mutate(
  
  enter = sample(lifetimes.empirical$enter, nrow(pred.df)),        # enter times distributed by our cap times
  sex = sample(lifetimes.empirical$sex, nrow(pred.df)),            # similar sex ratio
  mas = sample(lifetimes.empirical$mass, nrow(pred.df)),           # similar mass distribution
  hfl = 0,                                                         # all hfl at mean
  ret = 0,                                                         # control
  pil = 0)                                                         # control

# and construct the basis functions for the spline on hazard
# let's create a prediction df first
weeks.pred <- seq(1, 52)                 # weekly predictions

basis.pred <- bs(weeks.pred,       
                 knots = knot.list[-c(1, n.knots)],          
                 degree = 3,                     
                 intercept = FALSE)

# draw from posterior samples
sample.draws <- model.draws[sample(x = nrow(model.draws),
                            size = n.sim,
                            replace = TRUE), ]

# nested loop (draw, observation, follow-up)
lifetimes.all <- data.frame()

for (i in 1:n.sim) {
  
  # focal draw (drawn with replacement)
  focal.draw <- sample.draws %>% slice(i)
  
  # inner loop, by observation
  for (j in 1:nrow(pred.df)) {
    
    # subset
    focal.df <- pred.df %>% slice(j)
    
    # expand df to an entire year
    focal.df.1 <- data.frame(row = 1:n.steps)
    
    # fill other columns
    focal.df.1 <- focal.df.1 %>% 
      
      mutate(id = focal.df$id,
             week = recycle.time(focal.df$enter, n.steps),
             sex = focal.df$sex,
             mas = focal.df$mas,
             hfl = focal.df$hfl,
             ret = focal.df$ret,
             pil = focal.df$pil)
    
    # calculate predicted lambda
    # h0(t) (baseline hazard - spline)
    # spline weights
    w0 <- t(as.matrix(as.numeric(focal.draw[ , 2:7])))
    
    # multiply with 'sweep'
    w0.by.b <- sweep(basis.pred, 2, w0, `*`)
    
    # sum to create "normal" scale spline prediction
    w0.by.b.sum <- apply(w0.by.b, 1, sum)
    
    # join the baseline hazard prediction
    focal.df.2 <- focal.df.1 %>% 
      
      left_join(data.frame(week = 1:nrow(basis.pred),
                           h0 = exp(focal.draw$a0 + w0.by.b.sum)),
                by = join_by(week)) %>%
      
      # and calculate the lambda (total hazard)
      mutate(lambda = h0 * exp(focal.draw$b_sex * sex +
                               focal.draw$b_mas * mas +
                               focal.draw$b_hfl * hfl +
                               focal.draw$b_ret * ret +
                               focal.draw$b_pil * pil)) %>%
      
      # calculate the survival probability
      mutate(s = exp(-lambda)) %>%
      
      # and draw from Bernoulli trials using survival and censoring probs
      mutate(y = rbinom(nrow(focal.df.1), 1, p = 1 - s),
             y.cens = rbinom(nrow(focal.df.1), 1, p = 1 - censor.prob))
    
    # determine which comes first (if either) - death or censoring
    # case: neither event occurred
    if (sum(focal.df.2$y) == 0 & sum(focal.df.2$y.cens) == 0) {
      
      focal.lifetime <- nrow(focal.df.2)
      focal.died <- 0
      
    } 
    
    # case: death occurred, censoring did not
    if (sum(focal.df.2$y) > 0 & sum(focal.df.2$y.cens) == 0) {
      
      focal.lifetime <- which(focal.df.2$y == 1)[1]
      focal.died <- 1
      
    } 
    
    # case: death did not occur, censoring did 
    if (sum(focal.df.2$y) == 0 & sum(focal.df.2$y.cens) > 0) {
      
      focal.lifetime <- which(focal.df.2$y.cens == 1)[1]
      focal.died <- 0
      
    } 
    
    # case: death and censoring occurred
    if (sum(focal.df.2$y) > 0 & sum(focal.df.2$y.cens) > 0) {
      
      # did death happen first?
      which.first <- which(focal.df.2$y == 1)[1] < which(focal.df.2$y.cens == 1)[1]
      
      focal.lifetime <- ifelse(isTRUE(which.first),
                               which(focal.df.2$y == 1)[1],
                               which(focal.df.2$y.cens == 1)[1])
      
      focal.died <- ifelse(isTRUE(which.first),
                           1,
                           0)
      
    } 
    
    # bind into lifetime df
    focal.lifetime <- data.frame(draw = i,
                                 id = focal.df$id,
                                 lifetime = focal.lifetime,
                                 died = focal.died)
    
    lifetimes.all <- rbind(lifetimes.all, focal.lifetime)
    
  }
  
}

#_______________________________________________________________________________________________
# 4c. Examine empirical vs. simulated lifetime distributions ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # simulated distributions
  geom_density(data = lifetimes.all %>% filter(died == 1),
               aes(x = lifetime,
                   group = draw),
               color = "#33CCCC",
               fill = NA) +
  
  # empirical distribution
  geom_density(data = lifetimes.empirical %>% filter(died == 1),
               aes(x = lifetime),
               color = "black",
               fill = NA)
