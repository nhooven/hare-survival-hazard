# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Censor-as-dead scenarios
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Nov 2025 
# Date completed: 12 Nov 2025 
# Date last modified: 03 Nov 2026 
# R version: 4.5.2

# based on reviewer comments, this procedure is not well supported
# instead, let's fit a model (1 - base) where we don't flip any outcomes and
# another (2 - censor-as-dead) where we treat all censors as mortalities

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

fates <- read.csv("data_cleaned/fates_final_cleaned.csv")

#_______________________________________________________________________________
# 3. Clean ----
#_______________________________________________________________________________

fates.1 <- fates |>
  
  mutate(y.mort1 = ifelse(Event.type == "Mortality" & event == 1, 1, 0),
         y.mort2 = ifelse(Event.type == "Mortality" & event == 1 |
                          Event.type == "Censor" & General.cause == "Unknown" & event == 1, 1, 0)
         
  ) |>
  
  # keep only columns we need
  dplyr::select(
    
    deployment.1,
    cluster,
    site,
    MRID,
    Event.type,
    General.cause,
    y.mort1,
    y.mort2,
    week,
    study.week,
    Sex.1,
    Collar.type.1,
    p.dm,
    p.o,
    p.js,
    post1,
    post2,
    ret,
    pil
    
  ) %>%
  
  # rename as needed
  rename(sex = Sex.1,
         collar = Collar.type.1)

#_______________________________________________________________________________________________
# 9. Write to .csv ----
#_______________________________________________________________________________________________

write.csv(fates.1, "data_cleaned/fates_forModel.csv", row.names = F)
