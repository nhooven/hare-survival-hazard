# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08a - Posterior predictive checks - Bayesian p-value
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 27 Feb 2026 
# R version: 4.4.3

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

library(dplyr)       # manipulate and clean data
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
         study.week.s = (study.week - mean(study.week)) / sd(study.week),
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

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
  # same study.week, same sex, same cluster
  
  # function to take model draws and calculate the collective BLH
  calc_spline_pred <- function (x, y) {
    
    # y is the focal iteration
    # extract spline weights - wsc[s, c, b]
    w <- t(
      
      as.matrix(
        
        as.numeric(
          
          c(
            
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 1]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 2]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 3]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 4]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 5]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 6]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 7]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 8]")],
            y[paste0("wsc", "[", x$sex[1] + 1, ", ", x$cluster[1], ", 9]")]
            
            )
          
          )
        
        )
      
      )
    
    # multiply and sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
    
    # add to intercept and exponentiate
    blh = as.numeric(exp(y[paste0("a0sc[", x$sex[1] + 1, ", ", x$cluster[1], "]")] + w.by.b.sum[x$week[1]]))
    
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
                  log(y$hr_dm) * x$p.dm.s +
                  log(y$hr_open) * x$p.open.s +
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
    # must take the correct number of draws so our discrepancy doesn't explode!
    mutate(
      
      pois.draw.1 = ifelse(rpois(n(), haz.1) > 0, 1, 0),
      pois.draw.2 = ifelse(rpois(n(), haz.2) > 0, 1, 0),
      pois.draw.3 = ifelse(rpois(n(), haz.3) > 0, 1, 0)
      
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
#_______________________________________________________________________________________________

# loop through iterations
discrep.sum.ALL <- data.frame()
discrep.i.ALL <- data.frame()
  
for (i in 1:3000) {
  
  # extract focal iteration
  iter.draw1 <- model.fit.1[i, ]
  iter.draw2 <- model.fit.2[i, ]
  iter.draw3 <- model.fit.3[i, ]
  
  # loop through study.weeks
  discrep.i.j <- data.frame()
  
  for (j in 1:max(fates.1$study.week)) {
    
    # subset one week
    focal.week <- fates.1[fates.1$study.week == j, ]
    
    # apply baseline hazard calculation function
    # split by sex-cluster
    # create new variable
    focal.week$sc <- paste0(focal.week$sex + 1, "_", focal.week$cluster)
    
    sc.split <- split(focal.week, focal.week$sc)
    
    # returns a df including blh
    focal.week.1 <- do.call(rbind, lapply(sc.split, calc_blh))
    
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
    discrep.i.j <- rbind(discrep.i.j, discrep.df)
    
  } # j
  
  # sum over all study-weeks
  discrep.i.j.sum <- discrep.i.j %>%
    
    summarize(
      
      discrep.sim.1.sum = sum(discrep.sim.1),
      discrep.sim.2.sum = sum(discrep.sim.2),
      discrep.sim.3.sum = sum(discrep.sim.3),
      discrep.obs.1.sum = sum(discrep.obs.1),
      discrep.obs.2.sum = sum(discrep.obs.2),
      discrep.obs.3.sum = sum(discrep.obs.3)
      
      ) %>%
    
    # keep only the sum columns
    dplyr::select(discrep.sim.1.sum, discrep.sim.2.sum, discrep.sim.3.sum,
                  discrep.obs.1.sum, discrep.obs.2.sum, discrep.obs.3.sum) %>%
    
    # and add an index for i
    mutate(i = i)
    
  # boolean Dsim > Dobs
  discrep.i <- discrep.i.j.sum %>%
  
    mutate(
      
      Dsim.Dobs.1 = discrep.sim.1.sum > discrep.obs.1.sum,
      Dsim.Dobs.2 = discrep.sim.2.sum > discrep.obs.2.sum,
      Dsim.Dobs.3 = discrep.sim.3.sum > discrep.obs.3.sum
      
    ) %>%
    
    # keep only the Dsim > Dobs columns
    dplyr::select(Dsim.Dobs.1, Dsim.Dobs.2, Dsim.Dobs.3) %>%
    
    # and add an index for i
    mutate(i = i)
  
  # bind in both
  discrep.sum.ALL <- rbind(discrep.sum.ALL, discrep.i.j.sum)
  discrep.i.ALL <- rbind(discrep.i.ALL, discrep.i)
  
} # i
  

    
discrep.sum.ALL  # summed discrepancy values (for plotting)
discrep.i.ALL     # for Bayesian p-value


#_______________________________________________________________________________________________
# 6. Calculate Bayesian p-values ----
#_______________________________________________________________________________________________

(bayes.p <- discrep.all.i %>%
  
  summarize(scen1 = sum(Dsim.Dobs.1) / n(),
            scen2 = sum(Dsim.Dobs.2) / n(),
            scen3 = sum(Dsim.Dobs.3) / n())

)
