# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01a - Merge fates & covariate data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2023
# Date completed: 27 Nov 2023
# Date last modified: 09 Feb 2024 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(lubridate)       # work with dates
library(survival)        # survsplit function

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Raw data/fates_02_09_2024.csv")
covs <- read.csv("Raw data/covariates_02_09_2024.csv")

# define cutoff date
cutoff <- as.Date("2024-01-31", tz = "America/Los_Angeles")

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
                General.cause)

# covariates
covs.1 <- covs %>%
  
  # keep relevant columns
  dplyr::select(Site.ID,
                Animal.ID,
                Date,
                Ear.tag..,
                Collar.type,
                Sex,
                Mass,
                Bag.mass,
                HFL) %>%
  
  # rename columns
  rename(Site = Site.ID,
         Ear.tag = Ear.tag..)

#_______________________________________________________________________________________________
# 4. Format dates correctly ----
#_______________________________________________________________________________________________

# fates
fates.1 <- fates.1 %>%
  
  # format dates correctly
  mutate(Capture.date = as.Date(mdy(Capture.date)),
         Estimated.event.date = as.Date(mdy(Estimated.event.date))) %>%
  
  # change NAs in event date to 01-31-2024
  replace_na(list(Estimated.event.date = cutoff)) %>%
  
  # keep records with Capture.date < cutoff date
  filter(Capture.date < cutoff)

# covariates
covs.1 <- covs.1 %>% 
  
  # format dates correctly
  mutate(Date = as.Date(mdy(Date)))

#_______________________________________________________________________________________________
# 5. Split by week ----
#_______________________________________________________________________________________________

# here, we will use the survSplit function in the "survival" package to split the dataset by week
# FUTURE: Determine recurrent time horizon for biological year, not human year

# create numeric start and end variables
# define origin for our analysis (2022-01-01)
date.origin <- as.numeric(as.Date("2022-01-01"))

fates.1$start <- as.numeric(fates.1$Capture.date) - date.origin
fates.1$end <- as.numeric(fates.1$Estimated.event.date) - date.origin

# create vector of time points
cut.points <- seq(1,
                  as.numeric(cutoff) - date.origin,
                  7)

# create numeric status variable
fates.1$status.num <- ifelse(fates.1$Event.type == "Mortality",
                             1,
                             0)

# survSplit
fates.2 <- survSplit(Surv(start,
                          end,
                          status.num) ~ .,
                     data = fates.1,
                     cut = cut.points,
                     start = "start",
                     end = "end")

# keep relevant columns, create numeric week and year variables
fates.3 <- fates.2 %>%
  
  # select
  dplyr::select(Site,
                Animal.ID,
                Ear.tag,
                Collar.type,
                Sex,
                Event.type,
                start,
                end,
                status.num) %>%
  
  # create numeric week variable (back-transform)
  mutate(week = as.numeric(format(as.Date(end + date.origin, 
                                          origin = "1970-01-01"), 
                                  "%W"))) %>%
  
  # replace 0 weeks with 52, add year variable
  mutate(week = ifelse(week == 0, 
                       52, 
                       week),
         
         year = ifelse(end < 366,
                       2022,
                       ifelse(end > 730,
                              2024,
                              2023))) %>%
  
  # remove event type and start/end variables
  dplyr::select(-c(Event.type,
                   start,
                   end))

#_______________________________________________________________________________________________
# 6. Attribute covariate values to each individual/week observation ----
#_______________________________________________________________________________________________
# 6a. Calculate final mass ----
#_______________________________________________________________________________________________

covs.2 <- covs.1 %>%
  
  mutate(Final.mass = ifelse(is.na(Bag.mass),
                             Mass - 0.1,
                             Mass - Bag.mass)) %>%
  
  # remove Mass and Bag.mass
  dplyr::select(-c(Mass, 
                   Bag.mass))

#_______________________________________________________________________________________________
# 6b. Sex-stratified population means ----
#_______________________________________________________________________________________________

# determine mean mass and HFL
means <- covs.2 %>%
  
  # remove all individuals <= collaring mass
  filter(Final.mass >= 0.8) %>%
  
  # group by sex
  group_by(Sex) %>%
  
  # summarize
  summarize(mean.mass = round(mean(Final.mass), digits = 2),
            mean.hfl = round(mean(HFL, na.rm = TRUE), digits = 1))

#_______________________________________________________________________________________________
# 6c. Fill in missing values with "best-guess" values ----
#_______________________________________________________________________________________________

# NOTE: All observations MUST have an associated sex value

covs.3 <- data.frame()

for (i in unique(covs.2$Animal.ID)) {
  
  # subset covariate data
  focal.covs <- covs.2 %>% filter(Animal.ID == i) %>%
    
    # ensure that this is arranged by date
    arrange(Date) %>%
    
    # calculate a days column
    mutate(days = as.numeric(Date))
  
  # replace NAs
  # loop through all rows
  for (j in 1:nrow(focal.covs)) {
    
    # subset focal row
    focal.row <- focal.covs[j, ]
    
    # MASS
    # does the focal row have an NA for mass?
    if (is.na(focal.row$Final.mass)) {
      
      # are there any valid other values to use?
      # calculate day differences - temporary column that will be overwritten in each j
      focal.covs$days.diff <- abs(focal.row$days - focal.covs$days)
      
      if (any(!is.na(focal.covs$Final.mass) &   # must not be NA
              focal.covs$Final.mass >= 0.8 &    # mass must be above the cutoff
              focal.covs$days.diff > 0 &        # cannot be the focal row
              focal.covs$days.diff < 30)        # should be within 30 days
          ) {
        
        # use the value for the entry that minimizes the time difference
        focal.covs.temp <- focal.covs %>% filter(days.diff != 0 & is.na(Final.mass) == FALSE) 
        
        focal.row$Final.mass <- focal.covs.temp$Final.mass[focal.covs.temp$days.diff == min(focal.covs.temp$days.diff)]
      
        # if else, use the mean value by sex  
      } else {
        
        focal.row$Final.mass <- means$mean.mass[means$Sex == focal.row$Sex]
        
      }
      
      focal.row.1 <- focal.row %>% dplyr::select(Site,
                                                 Animal.ID,
                                                 Date,
                                                 Ear.tag,
                                                 Collar.type,
                                                 Sex,
                                                 HFL,
                                                 Final.mass,
                                                 days)
      
    } else {
      
      focal.row.1 <- focal.row %>% dplyr::select(Site,
                                                 Animal.ID,
                                                 Date,
                                                 Ear.tag,
                                                 Collar.type,
                                                 Sex,
                                                 HFL,
                                                 Final.mass,
                                                 days)
      
    }
    
    # HFL
    # does the focal row have an NA for HFL?
    if (is.na(focal.row$HFL)) {
      
      # are there any valid other values to use?
      # calculate day differences - temporary column that will be overwritten in each j
      focal.covs$days.diff <- abs(focal.row$days - focal.covs$days)
      
      if (any(!is.na(focal.covs$HFL) &          # must not be NA
              focal.covs$HFL >= 11.0 &          # HFL must be above the cutoff
              focal.covs$days.diff > 0 &        # cannot be the focal row
              focal.covs$days.diff < 30)        # should be within 30 days
      ) {
        
        # use the value for the entry that minimizes the time difference
        focal.covs.temp <- focal.covs %>% filter(days.diff != 0 & is.na(HFL) == FALSE)  
        
        focal.row$HFL <- focal.covs.temp$HFL[focal.covs.temp$days.diff == min(focal.covs.temp$days.diff)]
        
        # if else, use the mean value by sex  
      } else {
        
        focal.row$HFL <- means$mean.hfl[means$Sex == focal.row$Sex]
        
      }
      
      focal.row.1 <- focal.row %>% dplyr::select(Site,
                                                 Animal.ID,
                                                 Date,
                                                 Ear.tag,
                                                 Collar.type,
                                                 Sex,
                                                 HFL,
                                                 Final.mass,
                                                 days)
      
    } else {
      
      focal.row.1 <- focal.row %>% dplyr::select(Site,
                                                 Animal.ID,
                                                 Date,
                                                 Ear.tag,
                                                 Collar.type,
                                                 Sex,
                                                 HFL,
                                                 Final.mass,
                                                 days)
      
    }
    
    # bind into master df
    covs.3 <- rbind(covs.3, focal.row.1)
    
  }
  
}

# drop "days" column
covs.3 <- covs.3 %>% dplyr::select(-days)

#_______________________________________________________________________________________________
# 6d. Attribute covariates to all observations ----
#_______________________________________________________________________________________________

# remove any observations without Animal.IDs
fates.3 <- fates.3 %>% filter(Animal.ID != "")

# loop through all observations and attribute values from covs.1
fates.4 <- data.frame()

for (i in unique(fates.3$Animal.ID)) {
  
  # subset monitoring data
  focal.id <- fates.3 %>% filter(Animal.ID == i)
  
  # subset covariate data
  focal.covs <- covs.2 %>% filter(Animal.ID == i) %>%
    
    # ensure that this is arranged by date
    arrange(Date)
  
  # FOR NOW - assume the sex recorded in "fates" matches those in "covs"
  
  # fill in NAs
  # mass - any body mass < 0.8 kg should not have received a collar, so if a previous
  # handling occasion has a mass < than this cutoff, use a population mean (by sex)
  
  # nested loop to examine all Final.mass cells with NAs
  for (j in 1:nrow(focal.covs)) {
    
    focal.covs$Final.mass[j]
    
    # examine first entry
    if (focal.covs)
    
    if (focal.covs$Final.mass[j] >= 0.8) {
      
      
      
    }
    
  }
  
  
}