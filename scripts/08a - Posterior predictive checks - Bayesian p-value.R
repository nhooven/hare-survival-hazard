# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08a - Posterior predictive checks - Bayesian p-value
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 05 Nov 2026 
# R version: 4.5.2

#_______________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________

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

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(dplyr)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- readRDS("models/model_1.rds")

# dataset
fates <- read.csv("data_cleaned/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("data_cleaned/fates_final_cleaned_2.csv")

# BCI by deployment
fates.deploy <- readRDS("data_cleaned/fates_deploy.rds")
fates.deploy.df <- fates.deploy$fates.deploy

# day lookup table
day.lookup <- read.csv("data_cleaned/day_lookup.csv")

#_______________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________
# 3a. Convert mcmc.lists to dfs ----
#_______________________________________________________________________________

model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))

#_______________________________________________________________________________
# 3b. Attribute study-week to fates df ----
#_______________________________________________________________________________

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
  
  # join in BCI data
  left_join(fates.deploy.df %>% dplyr::select(deployment.1, BCI)) %>%
  
  # impute mean BCI 
  mutate(BCI = ifelse(is.na(BCI) == T,
                      fates.deploy$bci.mean,
                      BCI)) %>%
  
  # standardize variables
  mutate(BCI.s = (BCI - fates.deploy$bci.mean) / fates.deploy$bci.sd,
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________
# 3c. Spline for prediction ----
#_______________________________________________________________________________

# define knots
n.knots <- 9 + 1

# construct raw smooth (SAME as model fit)
smooth.raw <- smoothCon(s(week,
                          bs = "cp",        # cyclic penalized - what I really want
                          k = n.knots),
                        data = fates.1,
                        absorb.cons = T)    # identifiability constraint

# weeks to predict on
pred.df <- data.frame(week = 1:52)

# smooth2random
smooth.rand <- smooth2random(smooth.raw[[1]], 
                             vnames = "",
                             type = 2) # fixed and random

# new predict grid
smooth.pred <- PredictMat(smooth.raw[[1]], pred.df)

# new basis, correctly rotated
basis.pred <- smooth.pred %*% smooth.rand$trans.U %*% diag(smooth.rand$trans.D)

#_______________________________________________________________________________
# 4. Split dataset by study.week ----
#_______________________________________________________________________________

fates.1.split <- split(fates.1, fates.1$study.week)

#_______________________________________________________________________________
# 5. Calculate discrepancy ----
#_______________________________________________________________________________
# 5a. Helper functions ----
#_______________________________________________________________________________

# calculate baseline hazard
# we'll stratify individuals by their sex-cluster instead of calculating this separately
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
    # we only need one - these are all for the same week
    blh = exp(
      
      as.numeric(
        
        y[paste0("a0sc[", x$sex[1] + 1, ", ", x$cluster[1], "]")]
        
      ) + w.by.b.sum[x$week[1]]
      
    )
    
    # return
    return(blh)
    
  }
  
  # add to x
  x$blh <- calc_spline_pred(x, iter.draw1)
  
  # return x
  return(x)
  
}

# calculate full hazard
calc_haz <- function (x) {
  
  # x is each deployment from the focal.week
  # function to to take model draws and calculate full hazard
  calc_haz_1 <- function(x, y) {
    
    # y is the focal iteration
    # total hazard ratio prediction
    hr.pred = exp(y$b_bci * x$BCI.s +
                  y$b_col * x$collar +
                  y$b_dm * x$p.dm.s +
                  y$b_open * x$p.open.s +
                  y$b_post1 * x$post1 +
                  y$b_post2 * x$post2 +
                  y$b_ret_post1 * x$ret * x$post1 +
                  y$b_ret_post2 * x$ret * x$post2 +
                  y$b_pil_post1 * x$pil * x$post1 +
                  y$b_pil_post2 * x$pil * x$post2)
    
    # full hazard
    full.haz = x$blh * hr.pred
    
    # return
    return(full.haz)
    
  }
  
  x$haz <- calc_haz_1(x, iter.draw1)
  
  # return
  return(x)
  
}

# calculate expected mortalities and simulate
e_sim_morts <- function (x) {
  
  # x is the full focal week.2
  x.1 <- x %>%
    
    # keep only columns we need
    dplyr::select(y.mort1, haz) %>%
    
    # take Poisson draws using each hazard
    # must take the correct number of draws so our discrepancy doesn't explode!
    mutate(pois.draw = ifelse(rpois(n(), haz) > 0, 1, 0)) %>%
    
    # summarize
    summarize(
      
      # expected mort frequencies (sum the hazards)
      e.mort = sum(haz),
      
      # observed mort frequencies
      o.mort = sum(y.mort1),
      
      # simulated mort frequencies
      s.mort = sum(pois.draw)
      
    ) %>%
    
    # calculate the discrepancy measure
    # Dj = ((sim/obs events - predicted events)^2) / predicted events
    mutate(
      
      # simulated
      discrep.sim = ((s.mort - e.mort)^2) / e.mort,
      
      # observed
      discrep.obs = ((o.mort - e.mort)^2) / e.mort,
      
    ) %>%
    
    # keep only the discrep columns
    dplyr::select(discrep.sim, discrep.obs)
  
  # return
  return(x.1)
  
}

#_______________________________________________________________________________
# 5b. Loop ----
#_______________________________________________________________________________

# loop through iterations
discrep.sum.ALL <- data.frame()
discrep.i.ALL <- data.frame()
  
start.time <- Sys.time()

for (i in 1:3000) {
  
  # extract focal iteration
  iter.draw1 <- model.fit.1[i, ]
  
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
      
      discrep.sim.sum = sum(discrep.sim),
      discrep.obs.sum = sum(discrep.obs)
      
      ) %>%
    
    # keep only the sum columns
    dplyr::select(discrep.sim.sum, discrep.obs.sum) %>%
    
    # and add an index for i
    mutate(i = i)
    
  # boolean Dsim > Dobs
  discrep.i <- discrep.i.j.sum %>%
  
    mutate(
      
      Dsim.Dobs = discrep.sim.sum > discrep.obs.sum
      
    ) %>%
    
    # keep only the Dsim > Dobs columns
    dplyr::select(Dsim.Dobs) %>%
    
    # and add an index for i
    mutate(i = i)
  
  # bind in both
  discrep.sum.ALL <- rbind(discrep.sum.ALL, discrep.i.j.sum)
  discrep.i.ALL <- rbind(discrep.i.ALL, discrep.i)
  
  # timing
  if (i %% 50 == 0) {
    
    diff.time <- difftime(start.time, Sys.time(), units = "mins")
    
    print(
      
      paste0(
        
        "Completed iteration ", i, " of ", 3000, " - ", 
        round(as.numeric(diff.time), digits = 2), " mins"
        
      )
      
    )
    
  }
  
} # i

discrep.sum.ALL  # summed discrepancy values (for plotting)
discrep.i.ALL     # for Bayesian p-value

#_______________________________________________________________________________
# 6. Calculate Bayesian p-values ----
#_______________________________________________________________________________

(bayes.p <- discrep.i.ALL %>%
  
  summarize(mod1 = sum(Dsim.Dobs) / n())

)

# save
saveRDS(discrep.sum.ALL, "PPCs/discrep_sum_ALL.rds")
saveRDS(discrep.i.ALL, "PPCs/discrep_i_ALL.rds")
