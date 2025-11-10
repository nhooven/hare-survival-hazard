# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Merge fates & covariate data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2023
# Date completed: 27 Nov 2023
# Date last modified: 10 Nov 2025 
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
# 5. Dataset for fate classification submodel ----
#_______________________________________________________________________________________________

fate.class <- fates.1 %>%
  
  filter(is.na(Transmitter.lifetime) == FALSE)

# write to file
write.csv(fate.class, "Cleaned data/fate_class_cleaned.csv")

#_______________________________________________________________________________________________
# 6. Attribute covariate values to each individual/week observation ----
#_______________________________________________________________________________________________
# 6a. Read in covariate data and keep relevant columns ----
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
# 6b. Calculate final mass ----
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
# 6c. Fill in with matched measured values ----
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
fates.2.summary <- fates.2 %>%
  
  group_by(deployment)

sum(is.na(fates.2.summary$HFL))
sum(is.na(fates.2.summary$Final.mass))



# 11-10-2025
# Was going to impute covariates but ya know what, how about we still do this in the model?



#_______________________________________________________________________________________________
# 7. Standardize and format covariates ----
#_______________________________________________________________________________________________

fates.5 <- fates.4 %>%
  
  mutate(
    # use an indicator variable for sex (0 = F [intercept])
    Sex.1 = ifelse(Sex == "F",
                   0,
                   1),
    
    # use an indicator variable for Collar.type (0 = VHF-only)
    Collar.type.1 = ifelse(Collar.type == "VHF-ONLY",
                           0,
                           1),
    
    # center and scale continuous covariates
    Mass.1 = as.numeric(scale(Final.mass)),
    HFL.1 = as.numeric(scale(HFL)),
    
    # integer deployment
    deployment.1 = as.integer(as.factor(deployment))
  )








#_______________________________________________________________________________________________
# 5. Split by week and use correct identifiers for modeling ----
#_______________________________________________________________________________________________
# 5a. Define week cutoffs ----
#_______________________________________________________________________________________________

# here, we will use the survSplit function in the "survival" package to split the dataset by week

# define origin for our analysis (2022-10-01)
date.origin <- as.numeric(as.Date("2022-10-01"))

# create numeric start and end variables
fates.1 <- fates.1 %>%
  
  mutate(start = as.numeric(Capture.date) - date.origin,
         end = as.numeric(Estimated.event.date) - date.origin)

# create vector of weekly time points, from week 1 to the final week of the study
cut.points <- seq(1,
                  max(fates.1$end),
                  7)

#_______________________________________________________________________________________________
# 5b. Create lookup table for dates and weeks ----

# this table provides the date, numeric day, ordinal day of the year, integer month,
# calendar year, and integer week within the calendar year

#_______________________________________________________________________________________________

# final date of study
max.date = max(fates.1$Estimated.event.date)

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
# 5c. Create numeric status variables indicating end points ----

# this includes confirmed mortalities and informative censoring events

# we need to include dead transmitter events so we can include a collar dead vs. hare dead submodel

# this does not count:
# - collar removals
# - collar related morts after one week
# - slipped collars
# - any animals still alive at the end of the monitoring period 
# BUT we still want to include these observations

#_______________________________________________________________________________________________

fates.1 <- fates.1 %>% 
  
  mutate(status.num = ifelse(Event.type == "Mortality" |
                             (Event.type == "Censor" & 
                              General.cause %in% c("Unknown",
                                                   "Dead transmitter")),
                             1,
                             0))

#_______________________________________________________________________________________________
# 5d. Split dataset by cutoffs ----
#_______________________________________________________________________________________________

fates.2 <- survSplit(Surv(start,
                          end,
                          status.num) ~ .,
                     data = fates.1,
                     cut = cut.points,
                     start = "start",
                     end = "end")

# create numeric mort indicator variable
fates.2 <- fates.2 %>%
  
  mutate(mort = case_when(
    
    # weeks without an event get a zero
    Event.type == "Mortality" & status.num == 0 ~ 0,
    Event.type == "Censor" & status.num == 0 ~ 0,
    
    # weeks with a mortality get a 1
    Event.type == "Mortality" & status.num == 1 ~ 1,
    
    # weeks with an informative censoring event get an NA for imputation later
    Event.type == "Censor" & status.num == 1 ~ NA
    
  )
)

#_______________________________________________________________________________________________
# 5e. Keep relevant variables and create week variables ----
#_______________________________________________________________________________________________

# add correct "study.year.week" as a variable
fates.2$week <- NA

for (i in 1:length(fates.2$start)) {
  
  fates.2$week[i] <- day.lookup$study.year.week[which(day.lookup$study.day == fates.2$start[i])]
  
}

# clean
fates.3 <- fates.2 %>%
  
  # select columns we want
  dplyr::select(Site,
                Animal.ID,
                MRID,
                Capture.date,     # needed for splitting by deployment
                Collar.type,
                Sex,
                Event.type,
                General.cause,
                Specific.cause,
                start,
                end,
                mort,
                week) %>%
  
  # add year variable (from day.lookup)
  mutate(year = case_when(end < 93 ~ 2022,
                          end > 92 & end < 458 ~ 2023,
                          end > 457 & end < 824 ~ 2024,
                          end > 823 ~ 2025))

#_______________________________________________________________________________________________
# 6. Attribute covariate values to each individual/week observation ----
#_______________________________________________________________________________________________
# 6a. Read in covariate data and keep relevant columns ----
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
# 6b. Calculate final mass ----
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
# 6c. Fill in with matched measured values ----
#_______________________________________________________________________________________________

# split fates data by deployment
# new deployment column
fates.3$deployment = paste0(fates.3$MRID, "_", fates.3$Capture.date)

fates.3.split <- split(fates.3, f = fates.3$deployment)

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
fates.4 <- do.call(rbind, lapply(X = fates.3.split, FUN = attribute_covs))

# check how many missing observations we have - by deployment
fates.4.summary <- fates.4 %>%
  
  group_by(deployment) %>%
  
  slice(1)

sum(is.na(fates.4.summary$HFL))
sum(is.na(fates.4.summary$Final.mass))

#_______________________________________________________________________________________________
# 7. Standardize and format covariates ----
#_______________________________________________________________________________________________

fates.5 <- fates.4 %>%
  
  mutate(
         # use an indicator variable for sex (0 = F [intercept])
         Sex.1 = ifelse(Sex == "F",
                        0,
                        1),
         
         # use an indicator variable for Collar.type (0 = VHF-only)
         Collar.type.1 = ifelse(Collar.type == "VHF-ONLY",
                                0,
                                1),
         
         # center and scale continuous covariates
         Mass.1 = as.numeric(scale(Final.mass)),
         HFL.1 = as.numeric(scale(HFL)),
         
         # integer deployment
         deployment.1 = as.integer(as.factor(deployment))
         )

#_______________________________________________________________________________________________
# 8. Add treatment variable ----
#_______________________________________________________________________________________________

# new df to manipulate
fates.6 <- fates.5

# here, we will use three indicator variables to estimate the effect of each treatment
# compared to control (i.e., the intercept), while also accounting for pre- and post-treatment
# sensu Abele et al. (2013)

# week 1 for 2A, 2B, 3A, and 3B
# week 2 for 1A, 1B, 4A, and 4B

fates.6$post.trt <- NA
fates.6$trt.ret <- NA
fates.6$trt.pil <- NA

# assign 0 for pre-treatment and 1 for post-treatment
# 2 and 3
fates.6$post.trt[fates.6$Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
                 fates.6$week > 14 &
                 fates.6$year == 2023] <- 0

fates.6$post.trt[fates.6$Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
                 fates.6$week >= 1 &
                 fates.6$week <= 14 &
                 fates.6$year == 2023] <- 1

# 1 and 4
fates.6$post.trt[fates.6$Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
                 fates.6$week > 14 &
                 fates.6$week == 1 &
                 fates.6$year == 2023] <- 0

fates.6$post.trt[fates.6$Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
                 fates.6$week >= 2 &
                 fates.6$week <= 14 &
                 fates.6$year == 2023] <- 1

# all remaining NAs should be 1
fates.6$post.trt[is.na(fates.6$post.trt) == TRUE] <- 1

# assign treatments
# controls
fates.6$trt.ret[fates.6$Site %in% c("1C", "2C", "3C", "4C")] <- 0
fates.6$trt.pil[fates.6$Site %in% c("1C", "2C", "3C", "4C")] <- 0

# retention
fates.6$trt.ret[fates.6$Site %in% c("1A", "2B", "3B", "4A")] <- 1
fates.6$trt.pil[fates.6$Site %in% c("1A", "2B", "3B", "4A")] <- 0

# piling
fates.6$trt.ret[fates.6$Site %in% c("1B", "2A", "3A", "4B")] <- 0
fates.6$trt.pil[fates.6$Site %in% c("1B", "2A", "3A", "4B")] <- 1

#_______________________________________________________________________________________________
# 9. Add cluster variable (will be an index) ----
#_______________________________________________________________________________________________

fates.6$cluster <- as.numeric(substr(fates.6$Site, 1, 1))

#_______________________________________________________________________________________________
# 10. Write to csv ----
#_______________________________________________________________________________________________

write.csv(fates.6, "Cleaned data/fates_final_cleaned.csv")
