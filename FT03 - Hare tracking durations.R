# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Fate timing and incidence
# Script: FT03 - Hare tracking durations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Jan 2025
# Date completed: 
# Date last modified: 21 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(lubridate)       # work with dates

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_timing_cleaned_01_21_2025.csv")

# day lookup table (for reference)
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Monthly cutoffs ----
#_______________________________________________________________________________________________

month.cutoffs <- data.frame(breaks = c(1, 5.57, 9.86, 
                                       14.29, 18.71, 22.71, 
                                       27.14, 31.43, 35.86, 
                                       40.14, 44.57, 49.86),
                            labels = c("Oct", "Nov", "Dec", 
                                       "Jan", "Feb", "Mar", 
                                       "Apr", "May", "Jun", 
                                       "Jul", "Aug", "Sep"))

#_______________________________________________________________________________________________
# 4. Create and clean variables ----
#_______________________________________________________________________________________________
# 4a. Generic "event" variable ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  mutate(Event = case_when(Event.type == "Censor" & General.cause == "Removed" ~ "Removed",
                           Event.type == "Censor" & General.cause != "Removed" ~ "Censor",
                           Event.type == "On air" ~ "On air",
                           Event.type == "Mortality" ~ "Mortality"))
