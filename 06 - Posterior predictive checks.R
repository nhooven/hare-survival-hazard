# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06 - Posterior predictive checks
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 17 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# Bayesian p-value
# Here we'll apply a simplified procedure outlined in Jones et al. 2020
# For each week of our study, we'll use posterior draws to generate 
# predictions for our Poisson intensities, then calculate the frequency of 
# mortalities from simulated Poisson draws, then compare to observed frequencies

# we'll calculate the conditional (i.e., by week) difference between simulated and observed
# sum up all of these differences,
# do this for all iterations
# and look at how many iterations have more (conditional) simulated events than observed events
# should be close to 0.5

# Simulated lifetimes
# Here we'll simulate lifetimes from a subset of posterior draws, create Kaplan-Meier
# survivorship curves for each "population", then compare to the empirical K-M curve
# I'll just need to snap every individual to the same starting place

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- read.csv("Model outputs/model_1.csv")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________________________
# 3a. Attribute study-week to fates df ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  # join in year
  mutate(year = fates.year$year) %>%
  
  left_join(
    
    day.lookup %>%
      
      dplyr::select(
        
        year,
        study.week,
        study.year.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1))

#_______________________________________________________________________________________________
# 3b. Spline for prediction ----
#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# weeks to predict on
weeks.pred <- 1:52

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________________________
# 4. Split dataset by study.week ----
#_______________________________________________________________________________________________

fates.1.split <- split(fates.1, fates.1$study.week)

#_______________________________________________________________________________________________
# 5. Write Bayesian p-value function ----

# lapply() was a good thought, but I think just looping by iteration and doing the summation
# within the function is probably fine. We can also print the status by iteration

#_______________________________________________________________________________________________

bayes_p_mort <- function(
    
  draws = model.fit.1,
  response = "y.mort.scen1"
  
  ) {
  
  # functions
  calc_intensity <- function (y) {
    
    # BLH spline prediction
    # weights w
    w <- t(as.matrix(as.numeric(c(iter.draw[paste0("w", ".", y$sex_forest, "..1.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..2.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..3.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..4.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..5.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..6.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..7.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..8.")],
                                  iter.draw[paste0("w", ".", y$sex_forest, "..9.")]))))
    
    # multiply and sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
    
    # add to intercept and exponentiate
    blh = as.numeric(exp(iter.draw[paste0("a0.", y$sex_forest, ".")] + w.by.b.sum[y$week]))
    
    # coefficient prediction (must be on the log scale!)
    coef.pred = exp(log(iter.draw$hr_bci) * y$BCI.s +
                    log(iter.draw$hr_pil_total1) * y$post1 * y$pil +
                    log(iter.draw$hr_pil_total2) * y$post2 * y$pil +
                    log(iter.draw$hr_ret_total1) * y$post1 * y$ret +
                    log(iter.draw$hr_ret_total2) * y$post2 * y$ret)
    
    # total intensity
    intens = blh * coef.pred
    
    # Poisson draw
    pois.draw = rpois(n = 1, intens)
    
    # pack into a df
    y.intens.pois <- data.frame(intens = intens,
                                pois.draw = ifelse(pois.draw > 0, 1, pois.draw))
    
    return(y.intens.pois)
    
  }
  
  # loop through iterations
  discrep.sim.all.draws <- vector(length = nrow(draws))
  discrep.obs.all.draws <- vector(length = nrow(draws))
  
  for (i in 1:nrow(draws)) {
    
    # iteration
    iter.draw <- draws[i, ]
    
    # loop through study.weeks
    discrep.sim <- vector(length = max(fates.1$study.week))
    discrep.obs <- vector(length = max(fates.1$study.week))
    
    for (j in 1:max(fates.1$study.week)) {
      
      focal.week <- fates.1[fates.1$study.week == j, ]
      
      # function to apply through all observations
      indiv.split <- split(focal.week, focal.week$deployment.1)
      
      # apply function
      sim.events <- do.call(rbind, lapply(indiv.split, calc_intensity))
      
      # determine the expected number of events
      # this is just the sum of all intensities
      e.events = sum(sim.events$intens)
      
      # calculate the conditional discrepancy
      # Dj = ((sim/obs events - predicted events)^2) / predicted events
      discrep.sim[j] = sum(((sim.events$pois.draw - sim.events$intens)^2) / sim.events$intens)
      discrep.obs[j] = sum(((focal.week[[response]] - sim.events$intens)^2) / sim.events$intens)
      
    }
    
    # sum over all study.weeks
    discrep.sim.all.draws[i] = sum(discrep.sim)
    discrep.obs.all.draws[i] = sum(discrep.obs)
      
    # print status (every 100 iterations)
    
    if (i %% 5 == 0) {
      
      print(paste0("Completed iteration ", i, " of ", nrow(draws))) 
      
    }
    
  }
  
  # bind together
  discrep.sim.obs <- data.frame(sim = discrep.sim.all.draws,
                                obs = discrep.obs.all.draws)
  
  return(discrep.sim.obs)

}

#_______________________________________________________________________________________________
# 6. Apply function ----

# for all posterior iterations, this would take a long time. Best to leave it overnight

#_______________________________________________________________________________________________

bayes.diff <- bayes_p_mort(draws = draws, response = "y.mort.scen1")


bayes.diff <- data.frame(
  
  sim = discrep.sim.all.draws,
  obs = discrep.obs.all.draws,
  diff = discrep.sim.all.draws - discrep.obs.all.draws
  
  )

#_______________________________________________________________________________________________
# 7. Plot and calculate Bayes p-value ----
#_______________________________________________________________________________________________

length(which(bayes.diff$diff > 0)) / nrow(bayes.diff) # 0.4

# distribution
ggplot(data = bayes.diff,
       aes(x = diff)) +
  
  theme_bw() +
  
  geom_density(fill = "gray",
               alpha = 0.25) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  
  xlab("Dsim - Dobs")

# 1:1 plot
ggplot(data = bayes.diff,
       aes(x = obs,
           y = sim)) +
  
  theme_bw() +
  
  geom_point(shape = 21) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  xlab("Observed discrepancy") +
  ylab("Simulated discrepancy") +
  
  coord_cartesian(xlim = c(5000, 9000),
                  ylim = c(5000, 9000))

#_______________________________________________________________________________________________
# 8. Prepare dataset for Kaplan-Meier curves ----
#_______________________________________________________________________________________________
# 8a. Prepare empirical data ----

# Here we'll keep everything the same, but include a "t" variable denoting the week from the beginning of 
# monitoring per deployment

#_______________________________________________________________________________________________

# split by deployment
fates.1.deploy.split <- split(fates.1, fates.1$deployment.1)

# function to incorporate t since monitoring started
add_t <- function (x) {
  
  x$t = 1:nrow(x)
  
  return(x)
  
}

# apply
fates.2 <- do.call(rbind, lapply(fates.1.deploy.split, add_t))

#_______________________________________________________________________________________________
# 8b. Write function ----

# Ideally we can also use this within the simulation procedure
# we'll need to split data based on t since deployment

#_______________________________________________________________________________________________

kap_mei <- function (
    
  x,
  response = "y.mort.scen1"
  
  ) {
  
  # loop through all follow-up times
  all.S <- vector(length = max(x$t))
  s.i <- vector(length = max(x$t))
  
  all.S[1] = 1.0 # initialize at 100%
  s.i[1] = 1.0 # initialize at 100%
  
  for (i in 2:max(x$t)) {
    
    # subset
    indivs.i <- x[x$t == i, ]
    
    # calculate survival at interval
    s.i[i] = 1 - (sum(indivs.i[[response]]) / nrow(indivs.i))
    
    # cumulative survival at interval
    all.S[i] = prod(s.i[1:i])
    
  }
  
  # df for returning
  km.df <- data.frame(t = 1:max(x$t),
                      S = all.S)
  
  return(km.df)
  
}

km.df <- kap_mei(fates.2)

# plot
ggplot(data = km.df,
       aes(x = t,
           y = S)) +
  
  theme_bw() +
  
  geom_line() +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1))
