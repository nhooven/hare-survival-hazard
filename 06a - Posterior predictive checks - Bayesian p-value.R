# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06a - Posterior predictive checks - Bayesian p-value
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 01 Dec 2025 
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
model.fit.2 <- read.csv("Model outputs/model_2.csv")
model.fit.3 <- read.csv("Model outputs/model_3.csv")

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
        study.year.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1),
         study.week.s = (study.week - mean(study.week)) / sd(study.week))

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
# 5. Bayesian p-value ----

# lapply() was a good thought, but I think just looping by iteration and doing the summation
# within the function is probably fine. We can also print the status by iteration

# helper function 
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
                  log(iter.draw$hr_bci_study_week) * y$study.week.s * y$BCI.s +
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

#_______________________________________________________________________________________________

# this takes over an hour to run 1,000 iterations. Best to leave overnight

draws = model.fit.1

response = "y.mort.scen1"

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
    discrep.sim[j] = ((sum(sim.events$pois.draw) - sum(sim.events$intens))^2) / sum(sim.events$intens)
    discrep.obs[j] = ((sum(focal.week[[response]]) - sum(sim.events$intens))^2) / sum(sim.events$intens)
    
  }
  
  # sum over all study.weeks
  discrep.sim.all.draws[i] = sum(discrep.sim)
  discrep.obs.all.draws[i] = sum(discrep.obs)
    
  # print status (every 100 iterations)
  
  if (i %% 100 == 0) {
    
    print(paste0("Completed iteration ", i, " of ", nrow(draws))) 
    
  }
  
}

# bind together
bayes.diff <- data.frame(sim = discrep.sim.all.draws,
                         obs = discrep.obs.all.draws)

# 11-26-2025
# I stopped it here at 
nrow(bayes.diff[bayes.diff$sim > 0, ])
# 9183

# write to csv
write.csv(bayes.diff, "In progress data/bayes_diff_scen1.csv")

#_______________________________________________________________________________________________
# 6. calculate Bayes p-value ----
#_______________________________________________________________________________________________

bayes.diff$diff = bayes.diff$sim - bayes.diff$obs

length(which(bayes.diff$diff > 0)) / nrow(bayes.diff) # 0.12

#_______________________________________________________________________________________________
# 7. Plot ----
#_______________________________________________________________________________________________

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
  
  coord_cartesian(xlim = c(100, 250),
                  ylim = c(100, 250))

# slightly better
