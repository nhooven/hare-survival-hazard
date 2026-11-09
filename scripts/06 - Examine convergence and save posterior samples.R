# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06 - Examine convergence and save posterior samples
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 22 Jan 2026 
# Date completed: 22 Jan 2026 
# Date last modified: 09 Nov 2026 
# R version: 4.5.2

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(MCMCvis)
library(coda)
library(bayestestR)

#_______________________________________________________________________________
# 2. Read in posterior samples ----
#_______________________________________________________________________________

model.1 <- readRDS("models/model_clog.rds")

#_______________________________________________________________________________
# 3. Summarize main parameters ----

param.names <- colnames(model.1[[1]])

#_______________________________________________________________________________

MCMCsummary(model.1, params = c("a0s", "a0sc", param.names[c(13:22)]))

#_______________________________________________________________________________
# 4. Traceplots for betas ----
#_______________________________________________________________________________

# model 1
MCMCtrace(model.1, params = c(param.names[c(13:22)]), pdf = F)

#_______________________________________________________________________________
# 5. Save posterior summaries for betas ----
#_______________________________________________________________________________

model.1.betas <- MCMCsummary(model.1, params = c(param.names[c(13:22)]), hpd_prob = (0.9), HPD = T)

model.1.medians <- MCMCsummary(model.1, params = c(param.names[c(13:22)]), probs = 0.5)[, 3]

# write to tables
write.table(cbind(model.1.medians, model.1.betas), "clipboard", sep = "\t")

#_______________________________________________________________________________
# 6. Hazard ratios ----
#_______________________________________________________________________________

calc_hr <- function (x) {
  
  x.1 <- x %>%
    
    mutate(
      
      # easy ones
      hr.bci = exp(b_bci),
      hr.col = exp(b_col),
      hr.dm = exp(b_dm),
      hr.open = exp(b_open),
      
      # years alone
      hr.post1 = exp(b_post1),
      hr.post2 = exp(b_post2),
      
      # base treatment effects
      hr.ret.post1 = exp(b_ret_post1),
      hr.ret.post2 = exp(b_ret_post2),
      hr.pil.post1 = exp(b_pil_post1),
      hr.pil.post2 = exp(b_pil_post2),
      
      # total treatment effects
      hr.ret.post1.total = exp(b_post1 + b_ret_post1),
      hr.ret.post2.total = exp(b_post2 + b_ret_post2),
      hr.pil.post1.total = exp(b_post1 + b_pil_post1),
      hr.pil.post2.total = exp(b_post2 + b_pil_post2)
      
    ) %>%
  
    # pivot
    pivot_longer(hr.bci:hr.pil.post2.total) %>%
    
    # keep only relevant columns
    dplyr::select(name, value) %>%
    
    group_by(name) %>%
    
    summarize(
      
      med = median(value),
      sd = sd(value),
      lo.50 = as.numeric(hdi(value, ci = 0.50)[2]),
      up.50 = as.numeric(hdi(value, ci = 0.50)[3]),
      lo.90 = as.numeric(hdi(value, ci = 0.90)[2]),
      up.90 = as.numeric(hdi(value, ci = 0.90)[3])
      
    ) %>%
    
    ungroup()
  
  return(x.1)
  
}

model.fit.1 <- as.data.frame(do.call(rbind, model.1[ , c(13:22)]))

# use function
hr.1 <- calc_hr(model.fit.1)

# write.table
write.table(hr.1, "clipboard", sep = "\t")
