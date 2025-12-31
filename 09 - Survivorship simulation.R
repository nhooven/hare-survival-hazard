# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 09 - Survivorship simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 31 Dec 2025 
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# We'll create a sample of simulated individuals
  # even sex ratio
  # average (with variability) body condition
  # different exposure to treatment

# we could start three treatments at once OR
# wait until some midway period to implement treatment

# at least to start with, we'll only use a subsample of posterior draws for efficiency
# we'll use model 2 because it's a nice balance

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)
library(tictoc)          # timing

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.2 <- read.csv("Model outputs/model_2.csv")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Body condition score distribution ----
#_______________________________________________________________________________________________

fates.byIndiv <- fates %>%
  
  group_by(deployment.1) %>%
  
  slice(1)

hist(fates.byIndiv$BCI.1)

possible.BCI <- (fates.byIndiv$BCI.1 - mean(fates$BCI.1)) / sd(fates$BCI.1)

hist(possible.BCI)

#_______________________________________________________________________________________________
# 4. Simulated individual datasets ----

# each forest type / treatment combination will include 100 animals
# we'll simulate for 2 years (104 weeks)

#_______________________________________________________________________________________________
# 4a. Full exposure ----

# treatments here start at week 1

#_______________________________________________________________________________________________

# SET 1a - SFL
set.1a <- data.frame(
  
  indiv = 1:300,
  sex = rep(c(0, 1), each = 50, times = 3),
  BCI.1 = sample(possible.BCI, size = 300)
  
)

# sex_forest
set.1a$sex.forest = ifelse(set.1a$sex == 0, 1, 2)

# replicate rows
set.1a <- set.1a %>% slice(rep(1:n(), each = 104))
