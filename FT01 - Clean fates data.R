# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Fate timing and incidence
# Script: FT01 - Clean fates data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Jan 2025
# Date completed: 21 Jan 2025
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

fates <- read.csv("Raw data/fates_01_02_2025.csv")

# define cutoff date
cutoff <- as.Date("2024-12-31", tz = "America/Los_Angeles")
cutoff.char <- "12/31/2024"

# day lookup
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

# ensure that study.date is a date
day.lookup$study.date <- as.Date(ymd(day.lookup$study.date))

#_______________________________________________________________________________________________
# 3. Keep relevant columns ----
#_______________________________________________________________________________________________

# fates
fates.1 <- fates %>%
  
  # keep only records with include == "Y"
  filter(Include == "Y") %>%
  
  # keep relevant columns
  dplyr::select(Site,
                Animal.ID,
                Ear.tag,
                Collar.type,
                Sex,
                Capture.date,
                Estimated.event.date,
                Event.type,
                General.cause, 
                Specific.cause) %>%
  
  # if no event type (i.e., hare still on air at cutoff), add it
  mutate(Event.type = case_when(Event.type != "" ~ Event.type,
                                Event.type == "" ~ "On air"),
         Estimated.event.date = case_when(Estimated.event.date != "" ~ Estimated.event.date,
                                          Estimated.event.date == "" ~ cutoff.char))

#_______________________________________________________________________________________________
# 4. Format dates correctly ----
#_______________________________________________________________________________________________

# fates
fates.1 <- fates.1 %>%
  
  # format dates correctly
  mutate(Capture.date = as.Date(mdy(Capture.date)),
         Estimated.event.date = as.Date(mdy(Estimated.event.date))) %>%
  
  # change NAs in event date to cutoff
  replace_na(list(Estimated.event.date = cutoff)) %>%
  
  # keep records with Capture.date < cutoff date
  filter(Capture.date < cutoff)

#_______________________________________________________________________________________________
# 5. Add "week" variable ----
#_______________________________________________________________________________________________

# merge week into fates df
day.lookup.syw <- day.lookup %>% 
  
  dplyr::select(study.date,            # for merging with event date
                study.week,            # for defining treatment status
                study.year.week) %>%   # for plotting
  
  rename(Estimated.event.date = study.date)

fates.2 <- merge(fates.1, day.lookup.syw)

#_______________________________________________________________________________________________
# 6. Add treatment variables ----
#_______________________________________________________________________________________________

# new df to manipulate
fates.3 <- fates.2

# study.week 53 for 2A, 2B, 3A, and 3B
# study.week 54 for 1A, 1B, 4A, and 4B

# assign 0 for pre-treatment and 1 for post-treatment
fates.3 <- fates.3 %>%
  
  # assign 0 for pre-treatment and 1 for post-treatment
  mutate(post.trt = case_when(Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
                              study.week < 53 ~ 0,
                              Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
                              study.week >= 53 ~ 1,
                              Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
                              study.week < 54 ~ 0,
                              Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
                              study.week >= 54 ~ 1),
         # add which treatment
         trt = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ "control",
                         Site %in% c("1A", "2B", "3B", "4A") ~ "retention",
                         Site %in% c("1B", "2A", "3A", "4B") ~ "piling"))

#_______________________________________________________________________________________________
# 7. Add cluster variable (will be an index) ----
#_______________________________________________________________________________________________

fates.3$cluster <- as.numeric(substr(fates.3$Site, 1, 1))

#_______________________________________________________________________________________________
# 6. Write to .csv ----
#_______________________________________________________________________________________________

write.csv(fates.3, "Cleaned data/fates_timing_cleaned_01_21_2025.csv")
