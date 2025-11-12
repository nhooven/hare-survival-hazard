# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Merge fates & covariate data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2023
# Date completed: 27 Nov 2023
# Date last modified: 12 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(lubridate)       # work with dates
library(survival)        # survsplit function

#_______________________________________________________________________________________________
# 2. Read in fates data ----
#_______________________________________________________________________________________________

fates <- read.csv("Raw data/fates_final.csv")

#_______________________________________________________________________________________________
# 3. Keep relevant columns (fates) ----
#_______________________________________________________________________________________________

# fates
fates.1 <- fates %>%
  
  # keep only records with include == "Y"
  filter(Include == "Y") %>%
  
  # keep relevant columns
  dplyr::select(Site,
                Animal.ID,
                ET1,
                MRID,
                Collar.type,
                Sex,
                Capture.date,
                Estimated.event.date,
                Event.type,
                General.cause, 
                Specific.cause,
                Transmitter.lifetime)

#_______________________________________________________________________________________________
# 4. Format dates correctly ----
#_______________________________________________________________________________________________

# fates
fates.1 <- fates.1 %>%
  
  # format dates correctly
  mutate(Capture.date = as.Date(mdy(Capture.date)),
         Estimated.event.date = as.Date(mdy(Estimated.event.date)))

#_______________________________________________________________________________________________
# 5. Attribute covariate values to each individual/week observation ----
#_______________________________________________________________________________________________
# 5a. Read in covariate data and keep relevant columns ----
#_______________________________________________________________________________________________

covs <- read.csv("Raw data/covariates_final.csv")

covs.1 <- covs %>%
  
  # keep relevant columns
  dplyr::select(Site,
                AnimalID,
                Date,
                ET1,
                MRID,
                Collar.type,
                Sex,
                Total.mass..kg.,
                Bag.mass,
                Hind.foot.length..cm.,
                Weighed.w..collar.) %>%
  
  # rename columns
  rename(Mass = Total.mass..kg.,
         HFL = Hind.foot.length..cm.,
         WWC = Weighed.w..collar.) %>%
  
  # coerce date
  mutate(Date = as.Date(mdy(Date), tz = "America/Los_Angeles"))

#_______________________________________________________________________________________________
# 5b. Calculate final mass ----
#_______________________________________________________________________________________________

covs.2 <- covs.1 %>%
  
  mutate(Final.mass = case_when(
    
    # if bag mass is recorded, subtract it from measured mass
    is.na(Bag.mass) == F ~ Mass - Bag.mass,
    
    # if bag mass is not recorded, subtract the mean measured bag mass
    is.na(Bag.mass) == T ~ Mass - mean(Bag.mass, na.rm = T)
    
  )
  
  ) %>%
  
  # subtract the collar mass if the animal was weighed with a collar on
  # assume 40 g for all collars
  mutate(Final.mass = ifelse(
    
    WWC == "Y",
    Final.mass - 0.04,
    Final.mass
    
  )
  
  ) %>%
  
  # remove variables
  dplyr::select(-c(Mass, 
                   Bag.mass,
                   WWC))

#_______________________________________________________________________________________________
# 5c. Fill in with matched measured values ----
#_______________________________________________________________________________________________

# split fates data by deployment
# new deployment column
fates.1$deployment = paste0(fates.1$MRID, "_", fates.1$Capture.date)

fates.1.split <- split(fates.1, f = fates.1$deployment)

# function to take each deployment, find any corresponding capture entries,
# and attribute relevant mass + HFL observations
attribute_covs <- function (x) {
  
  # subset covs df
  indiv.covs <- covs.2 %>% filter(MRID == x$MRID[1])
  
  # subset all capture entries within 7 days (+/-) of the collaring event
  indiv.covs.cap.window <- indiv.covs %>% 
    
    filter(Date >= x$Capture.date[1] - 7 &
           Date <= x$Capture.date[1] + 7)
  
  # if such entries exist, use the measurements closest to the collaring event
  if (nrow(indiv.covs.cap.window) > 0) {
    
    # calculate days from collaring event
    indiv.covs.cap.window <- indiv.covs.cap.window %>%
      
      mutate(days.from.collaring = abs(x$Capture.date[1] - Date)) %>%
      
      # arrange by days from collaring
      arrange(days.from.collaring)
    
    # rank each set of measurements based upon days from collaring event, and
    # choose the first one that isn't NA (assuming that at least one isn't)
    # HFL
    if (sum(is.na(indiv.covs.cap.window$HFL)) < nrow(indiv.covs.cap.window)) {
      
      x$HFL = indiv.covs.cap.window$HFL[which.min(is.na(indiv.covs.cap.window$HFL))]
      
    } else {
      
      x$HFL = NA
      
    }
    
    # Final.mass
    if (sum(is.na(indiv.covs.cap.window$Final.mass)) < nrow(indiv.covs.cap.window)) {
      
      x$Final.mass = indiv.covs.cap.window$Final.mass[which.min(is.na(indiv.covs.cap.window$Final.mass))]
      
    } else {
      
      x$Final.mass = NA
      
    }
    
    
    # if such entries do not exist, assign NAs
  } else {
    
    x$HFL = NA
    x$Final.mass = NA
    
  }
  
  return(x)
  
}

# apply function
fates.2 <- do.call(rbind, lapply(X = fates.1.split, FUN = attribute_covs))

# check how many missing observations we have - by deployment
sum(is.na(fates.2$HFL))
sum(is.na(fates.2$Final.mass))

# format covariates
fates.3 <- fates.2 %>%
  
  mutate(
    # use an indicator variable for sex (0 = F [intercept])
    Sex.1 = ifelse(Sex == "F",
                   0,
                   1),
    
    # use an indicator variable for Collar.type (0 = VHF-only)
    Collar.type.1 = ifelse(Collar.type == "VHF-ONLY",
                           0,
                           1),
    
    # integer deployment
    deployment.1 = as.integer(as.factor(deployment))
  )

#_______________________________________________________________________________________________
# 6. Split by week and use correct identifiers for modeling ----
#_______________________________________________________________________________________________
# 6a. Define week cutoffs ----
#_______________________________________________________________________________________________

# here, we will use the survSplit function in the "survival" package to split the dataset by week

# define origin for our analysis (2022-10-01)
date.origin <- as.numeric(as.Date("2022-10-01"))

# create numeric start and end variables
fates.3 <- fates.3 %>%
  
  mutate(start = as.numeric(Capture.date) - date.origin,
         end = as.numeric(Estimated.event.date) - date.origin)

# create vector of weekly time points, from week 1 to the final week of the study
cut.points <- seq(1,
                  max(fates.3$end),
                  7)

#_______________________________________________________________________________________________
# 6b. Create lookup table for dates and weeks ----

# this table provides the date, numeric day, ordinal day of the year, integer month,
# calendar year, and integer week within the calendar year

#_______________________________________________________________________________________________

# final date of study
max.date = max(fates.3$Estimated.event.date)

# create lookup table
day.lookup <- data.frame(study.date = seq(as.Date("2022-10-01"), max.date, by = 1),
                         study.day = 1:length(seq(as.Date("2022-10-01"), max.date, by = 1)),
                         year.day = yday(seq(as.Date("2022-10-01"), max.date, by = 1)),
                         month = month(seq(as.Date("2022-10-01"), max.date, by = 1)),
                         year = year(seq(as.Date("2022-10-01"), max.date, by = 1)),
                         year.week = week(seq(as.Date("2022-10-01"), max.date, by = 1)))

# add in study week identifier
week.id <- rep(1:length(cut.points), each = 7)

day.lookup$study.week <- week.id[1:(length(week.id) - 
                                   (length(week.id) - nrow(day.lookup)))]

# add study-year-week (what we'll use for modeling)
# this begins again at the start of Oct each year (i.e., the recurrent time origin/horizon)
study.year.week.id <- c(rep(1:52, each = 7), rep(1:52, each = 7),
                        rep(1:52, each = 7), rep(1:52, each = 7))

day.lookup$study.year.week <- study.year.week.id[1:(length(study.year.week.id) - 
                                                   (length(study.year.week.id) - nrow(day.lookup)))]

# write to .csv for use later
write.csv(day.lookup, "Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 6c. Create event variable ----

# Mortalities and unknown censor events get a 1
# all other observations get a zero
# this is so we can easily retrieve what we need to predict on

#_______________________________________________________________________________________________

fates.3 <- fates.3 %>% 
  
  mutate(
    
    event = ifelse(Event.type == "Mortality" |
                   (Event.type == "Censor" & 
                   General.cause == "Unknown"),
                   1,
                   0),
    
  )
                    
#_______________________________________________________________________________________________
# 6d. Split dataset by events ----
#_______________________________________________________________________________________________

fates.4 <- survSplit(Surv(start,
                          end,
                          event) ~ .,
                     data = fates.3,
                     cut = cut.points,
                     start = "start",
                     end = "end")

#_______________________________________________________________________________________________
# 6e. Mortality indicators ----

# SCENARIO 1 (all unknown censors remain zero) - y.mort.scen1
# ALL mortality events get a 1
# all other observations get a zero

# y.pred.scen1
# confirmed predation and unknown mortalities get a 1
# all other observations get a zero (including confirmed other mort types)

#_______________________________________________________________________________________________

fates.4 <- fates.4 %>% 
  
  mutate(
    
    y.mort.scen1 = ifelse(event == 1 & Event.type == "Mortality",
                          1,
                          0),
    
    y.pred.scen1 = ifelse(event == 1 & 
                          Event.type == "Mortality" & 
                          (General.cause == "Predation" | General.cause == "Unknown"),
                          1,
                          0)
    
  )

#_______________________________________________________________________________________________
# 6f. Create week variable and keep only columns we need ----
#_______________________________________________________________________________________________

# add correct "study.year.week" as a variable
fates.4$week <- NA

for (i in 1:length(fates.4$start)) {
  
  fates.4$week[i] <- day.lookup$study.year.week[which(day.lookup$study.day == fates.4$start[i])]
  
}

# clean
fates.5 <- fates.4 %>%
  
  # select columns we want
  dplyr::select(Site,
                Animal.ID,
                MRID,
                deployment.1,
                Event.type,
                General.cause,
                Specific.cause,
                start,
                end,
                event,
                y.mort.scen1,
                y.pred.scen1,
                week,
                Sex.1,
                Collar.type.1,
                HFL,
                Final.mass,
                Transmitter.lifetime) %>%
  
  # add year variable (from day.lookup)
  mutate(year = case_when(end < 93 ~ 2022,
                          end > 92 & end < 458 ~ 2023,
                          end > 457 & end < 824 ~ 2024,
                          end > 823 ~ 2025))

#_______________________________________________________________________________________________
# 7. Add treatment variable ----
#_______________________________________________________________________________________________

# new df to manipulate
fates.6 <- fates.5

# here, we will use three indicator variables to estimate the effect of each treatment
# compared to control (i.e., the intercept), while also accounting for  each year:
# pre, post1, and post2
# sensu Abele et al. (2013)

# week 1 for 2A, 2B, 3A, and 3B
# week 2 for 1A, 1B, 4A, and 4B

# assign indicators for each period
fates.7 <- fates.6 %>%
  
  # join in study-year-week to make this easy
  left_join(
    
    day.lookup %>% 
      
      dplyr::select(
        
        year,
        study.year.week,
        study.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1) %>%
      
      ungroup()
    
  ) %>%
  
  # use study.week to assign indicators
  mutate(
    
    # post 1
    post1 = case_when(
      
      # 2 and 3
      # pre-treatment
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week < 53 ~ 0,
      
      # post-treatment 1
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week >= 53 &
      study.week <= 105 ~ 1,
      
      # post-treatment 2
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week > 105 ~ 0,
      
      # 1 and 4
      # pre-treatment
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week <= 53 ~ 0,
      
      # post-treatment 1
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week >= 54 &
      study.week <= 106 ~ 1,
      
      # post-treatment 2
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week > 106 ~ 0
      
       ),
    
    post2 = case_when(
      
      # 2 and 3
      # pre-treatment
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week < 53 ~ 0,
      
      # post-treatment 1
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week >= 53 &
      study.week <= 105 ~ 0,
      
      # post-treatment 2
      Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
      study.week > 105 ~ 1,
      
      # 1 and 4
      # pre-treatment
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week <= 53 ~ 0,
      
      # post-treatment 1
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week >= 54 &
      study.week <= 106 ~ 0,
      
      # post-treatment 2
      Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
      study.week > 106 ~ 1
      
    ),
    
    # assign treatments
    ret = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ 0,
                    Site %in% c("1A", "2B", "3B", "4A") ~ 1,
                    Site %in% c("1B", "2A", "3A", "4B") ~ 0),
    
    pil = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ 0,
                    Site %in% c("1A", "2B", "3B", "4A") ~ 0,
                    Site %in% c("1B", "2A", "3A", "4B") ~ 1)
    
  ) %>%
  
  # remove study.week
  dplyr::select(-study.week)


#_______________________________________________________________________________________________
# 9. Add index cluster site variables ----
#_______________________________________________________________________________________________

fates.8 <- fates.7 %>% 
  
  mutate(cluster = as.numeric(substr(Site, 1, 1)),
         site = as.integer(factor(Site)))

#_______________________________________________________________________________________________
# 10. Write to csv ----
#_______________________________________________________________________________________________

write.csv(fates.8, "Cleaned data/fates_final_cleaned.csv")
