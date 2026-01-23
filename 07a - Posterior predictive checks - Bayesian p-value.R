# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 07a - Posterior predictive checks - Bayesian p-value
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 23 Jan 2026 
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
model.fit.1 <- readRDS("Model outputs/model_1.rds")
model.fit.2 <- readRDS("Model outputs/model_2.rds")
model.fit.3 <- readRDS("Model outputs/model_3.rds")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________________________
# 3a. Convert mcmc.lists to dfs ----
#_______________________________________________________________________________________________

model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

#_______________________________________________________________________________________________
# 3b. Attribute study-week to fates df ----
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
# 3c. Spline for prediction ----
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
# 5. Calculate discrepancy ----
#_______________________________________________________________________________________________
# 5a. Helper functions ----
#_______________________________________________________________________________________________

# calculate baseline hazard
# we'll stratify individuals by their sex-ft instead of calculating this separately
calc_blh <- function (x) {
  
  # x is a collection of individuals subject to the same BLH
  # same study.week, same sex, same forest type
  
  # function to take model draws and calculate the collective BLH
  calc_spline_pred <- function (x, y) {
    
    # y is the focal iteration
    # extract spline weights w
    w <- t(as.matrix(as.numeric(c(y[paste0("w", "[", x$sex_forest[1], ", 1]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 2]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 3]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 4]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 5]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 6]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 7]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 8]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 9]")]))))
    
    # multiply and sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
    
    # add to intercept and exponentiate
    blh = as.numeric(exp(y[paste0("a0[", x$sex_forest[1], "]")] + w.by.b.sum[x$week[1]]))
    
    # return
    return(blh)
    
  }
  
  # add to x
  x$blh.1 <- calc_spline_pred(x, iter.draw1)
  x$blh.2 <- calc_spline_pred(x, iter.draw2)
  x$blh.3 <- calc_spline_pred(x, iter.draw3)
  
  # return x
  return(x)
  
}

# calculate full hazard
calc_haz <- function (x) {
  
  # x is each deployment from the focal.week
  # function to to take model draws and calculate full hazard
  calc_haz_1 <- function(x, y, blh = "blh.1") {
    
    # y is the focal iteration
    # total hazard ratio prediction
    hr.pred = exp(log(y$hr_bci) * x$BCI.s +
                    log(y$hr_bci_study_week) * x$study.week.s * x$BCI.s +
                    log(y$hr_pil_total1) * x$post1 * x$pil +
                    log(y$hr_pil_total2) * x$post2 * x$pil +
                    log(y$hr_ret_total1) * x$post1 * x$ret +
                    log(y$hr_ret_total2) * x$post2 * x$ret)
    
    # full hazard
    full.haz = x[ , blh] * hr.pred
    
    # return
    return(full.haz)
    
  }
  
  x$haz.1 <- calc_haz_1(x, iter.draw1, "blh.1")
  x$haz.2 <- calc_haz_1(x, iter.draw2, "blh.2")
  x$haz.3 <- calc_haz_1(x, iter.draw3, "blh.3")
  
  # return
  return(x)
  
}

# calculate expected mortalities and simulate
e_sim_morts <- function (x) {
  
  # x is the full focal week.2
  x.1 <- x %>%
    
    # keep only columns we need
    dplyr::select(y.mort.scen1, y.mort.scen2, y.mort.scen3,
                  haz.1, haz.2, haz.3) %>%
    
    # take Poisson draws using each hazard
    mutate(
      
      pois.draw.1 = rpois(1, haz.1),
      pois.draw.2 = rpois(1, haz.2),
      pois.draw.3 = rpois(1, haz.3)
      
    ) %>%
    
    # summarize
    summarize(
      
      # expected mort frequencies (sum the hazards)
      e.mort.1 = sum(haz.1),
      e.mort.2 = sum(haz.2),
      e.mort.3 = sum(haz.3),
      
      # observed mort frequencies
      o.mort.1 = sum(y.mort.scen1),
      o.mort.2 = sum(y.mort.scen2),
      o.mort.3 = sum(y.mort.scen3),
      
      # simulated mort frequencies
      s.mort.1 = sum(pois.draw.1),
      s.mort.2 = sum(pois.draw.2),
      s.mort.3 = sum(pois.draw.3)
      
    ) %>%
    
    # calculate the discrepancy measure
    # Dj = ((sim/obs events - predicted events)^2) / predicted events
    mutate(
      
      # simulated
      discrep.sim.1 = ((s.mort.1 - e.mort.1)^2) / e.mort.1,
      discrep.sim.2 = ((s.mort.2 - e.mort.2)^2) / e.mort.2,
      discrep.sim.3 = ((s.mort.3 - e.mort.3)^2) / e.mort.3,
      
      # observed
      discrep.obs.1 = ((o.mort.1 - e.mort.1)^2) / e.mort.1,
      discrep.obs.2 = ((o.mort.2 - e.mort.2)^2) / e.mort.2,
      discrep.obs.3 = ((o.mort.3 - e.mort.3)^2) / e.mort.3
      
    ) %>%
    
    # keep only the discrep columns
    dplyr::select(discrep.sim.1:discrep.obs.3)
  
  # return
  return(x.1)
  
}

#_______________________________________________________________________________________________
# 5b. Loop ----

# one iteration takes 6.16 s
(6.16 * 3000) / 3600
# just over 5 hours

#_______________________________________________________________________________________________

# loop through iterations
discrep.all.i <- data.frame()

for (i in 1:3000) {
  
  # extract focal iteration
  iter.draw1 <- model.fit.1[i, ]
  iter.draw2 <- model.fit.2[i, ]
  iter.draw3 <- model.fit.3[i, ]
  
  # loop through study.weeks
  discrep.df.all.j <- data.frame()
  
  for (j in 1:max(fates.1$study.week)) {
    
    # subset one week
    focal.week <- fates.1[fates.1$study.week == j, ]
    
    # apply baseline hazard calculation function
    # split by sex-FT
    sft.split <- split(focal.week, focal.week$sex_forest)
    
    # returns a df including blh
    focal.week.1 <- do.call(rbind, lapply(sft.split, calc_blh))
    
    # apply full hazard calculation function
    # split by deployment
    deploy.split <- split(focal.week.1, focal.week.1$deployment.1)
    
    # returns a df with full hazard columns
    focal.week.2 <- do.call(rbind, lapply(deploy.split, calc_haz))
    
    # calculate expected mortalities and simulate
    # returns a df with scen and boolean column Dsim > Dobs
    # we'll add a j index
    discrep.df <- cbind(j, e_sim_morts(focal.week.2))
    
    # bind into master df
    discrep.df.all.j <- rbind(discrep.df.all.j, discrep.df)
    
  }
  
  # sum over all study-weeks
  discrep.j.sum <- discrep.df.all.j %>%
    
    summarize(discrep.sim.1.sum = sum(discrep.sim.1),
              discrep.sim.2.sum = sum(discrep.sim.2),
              discrep.sim.3.sum = sum(discrep.sim.3),
              discrep.obs.1.sum = sum(discrep.obs.1),
              discrep.obs.2.sum = sum(discrep.obs.2),
              discrep.obs.3.sum = sum(discrep.obs.3)) %>%
    
    # boolean Dsim > Dobs
    mutate(
      
      Dsim.Dobs.1 = discrep.sim.1.sum > discrep.obs.1.sum,
      Dsim.Dobs.2 = discrep.sim.2.sum > discrep.obs.2.sum,
      Dsim.Dobs.3 = discrep.sim.3.sum > discrep.obs.3.sum
      
    ) %>%
    
    # keep only the Dsim > Dobs columns
    dplyr::select(Dsim.Dobs.1, Dsim.Dobs.2, Dsim.Dobs.3) %>%
    
    # and add an index for i
    mutate(i = i)
  
  # and bind in
  discrep.all.i <- rbind(discrep.all.i, discrep.j.sum)
  
  # print status and save to .csv every 50 iterations
  if (i %% 50 == 0) {
    
    print(paste0("Completed iteration ", i, " of 3000"))
    
    write.csv(discrep.all.i, "discrep_all_i.csv")
    
  }
    
}

# 01-23-2026
# cleaned up a lot of code. Now to run it on the lab pc





#_______________________________________________________________________________________________
# 6. Calculate Bayes p-value ----
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
