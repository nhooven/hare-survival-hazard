# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06 - Examine convergence and save posterior samples
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 22 Jan 2026 
# Date completed: 22 Jan 2026 
# Date last modified: 12 Mar 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulateand clean data
library(MCMCvis)
library(coda)
library(bayestestR)

#_______________________________________________________________________________________________
# 2. Read in posterior samples ----
#_______________________________________________________________________________________________

model.1 <- readRDS("Model outputs/model_1.rds")
model.2 <- readRDS("Model outputs/model_2.rds")
model.3 <- readRDS("Model outputs/model_3.rds")

#_______________________________________________________________________________________________
# 3. Summarize main parameters ----

param.names <- colnames(model.1[[1]])

#_______________________________________________________________________________________________

# model 1
MCMCsummary(model.1, params = c("a0s", "a0sc", param.names[c(13:21)]))

# model 2
MCMCsummary(model.2, params = c("a0s", "a0sc", param.names[c(13:21)]))

# model 3
MCMCsummary(model.3, params = c("a0s", "a0sc", param.names[c(13:21)]))

#_______________________________________________________________________________________________
# 4. Traceplots for betas ----
#_______________________________________________________________________________________________

# model 1
MCMCtrace(model.1, params = c(param.names[c(13:21)]), pdf = F)

# model 2
MCMCtrace(model.2, params = c(param.names[c(13:21)]), pdf = F)

# model 3
MCMCtrace(model.3, params = c(param.names[c(13:21)]), pdf = F)

#_______________________________________________________________________________________________
# 5. Save posterior summaries for betas ----
#_______________________________________________________________________________________________

model.1.betas <- MCMCsummary(model.1, params = c(param.names[c(13:21)]), hpd_prob = (0.9), HPD = T)
model.2.betas <- MCMCsummary(model.2, params = c(param.names[c(13:21)]), hpd_prob = (0.9), HPD = T)
model.3.betas <- MCMCsummary(model.3, params = c(param.names[c(13:21)]), hpd_prob = (0.9), HPD = T)

# write to tables
write.table(model.1.betas, "clipboard", sep = "\t")
write.table(model.2.betas, "clipboard", sep = "\t")
write.table(model.3.betas, "clipboard", sep = "\t")

#_______________________________________________________________________________________________
# 6. Hazard ratios ----
#_______________________________________________________________________________________________

calc_hr <- function (x) {
  
  x.1 <- x %>%
    
    mutate(
      
      # easy ones
      hr.bci = exp(b_bci),
      hr.dm = exp(b_dm),
      hr.open = exp(b_open),
      
      # years alone
      hr.post1 = exp(b_post1),
      hr.post2 = exp(b_post2),
      
      # total treatment effects
      hr.ret.post1.total = exp(b_post1 + b_ret_post1),
      hr.ret.post2.total = exp(b_post2 + b_ret_post2),
      hr.pil.post1.total = exp(b_post1 + b_pil_post1),
      hr.pil.post2.total = exp(b_post2 + b_pil_post2),
      
      # year comparison effects
      hr.ret.post1.comp = exp(b_ret_post1 - b_post1),
      hr.ret.post2.comp = exp(b_ret_post2 - b_post2),
      hr.pil.post1.comp = exp(b_pil_post1 - b_post1),
      hr.pil.post2.comp = exp(b_pil_post2 - b_post2)
      
    ) %>%
  
    # pivot
    pivot_longer(hr.bci:hr.pil.post2.comp) %>%
    
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

model.fit.1 <- as.data.frame(do.call(rbind, model.1[ , c(13:21)]))
model.fit.2 <- as.data.frame(do.call(rbind, model.2[ , c(13:21)]))
model.fit.3 <- as.data.frame(do.call(rbind, model.3[ , c(13:21)]))

# use function
hr.1 <- calc_hr(model.fit.1)
hr.2 <- calc_hr(model.fit.2)
hr.3 <- calc_hr(model.fit.3)

# write.table
write.table(hr.1, "clipboard", sep = "\t")
write.table(hr.2, "clipboard", sep = "\t")
write.table(hr.3, "clipboard", sep = "\t")
