# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 05 - Examine convergence
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 22 Jan 2026 
# Date completed: 22 Jan 2026 
# Date last modified: 26 Jan 2026 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(MCMCvis)

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
MCMCsummary(model.1, params = c("a0", "lambda", param.names[c(7:12)]))

# model 2
MCMCsummary(model.2, params = c("a0", "lambda", param.names[c(7:12)]))

# model 3
MCMCsummary(model.3, params = c("a0", "lambda", param.names[c(7:12)]))

#_______________________________________________________________________________________________
# 4. Traceplots for HRs ----
#_______________________________________________________________________________________________

# model 1
MCMCtrace(model.1, params = c(param.names[c(7:12)]), pdf = F)

# model 2
MCMCtrace(model.2, params = c(param.names[c(7:12)]), pdf = F)

# model 3
MCMCtrace(model.3, params = c(param.names[c(7:12)]), pdf = F)

#_______________________________________________________________________________________________
# 5. Save Rhat and ESS for hazard ratios ----
#_______________________________________________________________________________________________

write.table(MCMCsummary(model.1, params = c(param.names[c(7:12)])), "clipboard", sep = "\t")
write.table(MCMCsummary(model.2, params = c(param.names[c(7:12)])), "clipboard", sep = "\t")
write.table(MCMCsummary(model.3, params = c(param.names[c(7:12)])), "clipboard", sep = "\t")
