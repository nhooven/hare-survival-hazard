# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Merge fates & covariate data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2023
# Date completed: 27 Nov 2023
# Date last modified: 30 Aug 2024 
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

fates <- read.csv("Raw data/fates_08_30_2024.csv")
covs <- read.csv("Raw data/covariates_08_30_2024.csv")

# define cutoff date
cutoff <- as.Date("2024-08-31", tz = "America/Los_Angeles")

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
  
  # change NAs in event date to cutoff
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

# create numeric status variable (indicating if this is an end point or not)
# here we should only indicate censoring events that are UNINFORMATIVE
fates.1 <- fates.1 %>% 
  
  mutate(status.num = ifelse(Event.type == "Mortality" |
                            (Event.type == "Censor" & 
                             General.cause %in% c("Unknown",
                                                  "Dead transmitter")),
                             1,
                             0))

# survSplit
fates.2 <- survSplit(Surv(start,
                          end,
                          status.num) ~ .,
                     data = fates.1,
                     cut = cut.points,
                     start = "start",
                     end = "end")

# change NAs in "start" to the stop time and create numeric indicator variables
fates.2 <- fates.2 %>%
  
  mutate(start = ifelse(is.na(start) == TRUE,
                        end,
                        start)) %>%
  
  mutate(mort = ifelse(Event.type == "Mortality" &
                       status.num == 1,
                       1,
                       0),
         cens = ifelse(Event.type == "Censor" &
                       status.num == 1,
                       1,
                       0))

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
                mort,
                cens) %>%
  
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
  dplyr::select(-c(Event.type))

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

# add numeric measurement day
covs.4 <- covs.3 %>% mutate(measure.date = as.numeric(Date) - date.origin)

#_______________________________________________________________________________________________
# 6d. Attribute covariates to all observations ----
#_______________________________________________________________________________________________

# remove any observations without Animal.IDs
fates.3 <- fates.3 %>% filter(Animal.ID != "")

# loop through all observations and attribute values from covs.3
fates.4 <- data.frame()

for (i in unique(fates.3$Animal.ID)) {
  
  # subset monitoring data
  focal.id <- fates.3 %>% filter(Animal.ID == i)
  
  # subset covariate data
  focal.covs <- covs.4 %>% filter(Animal.ID == i) %>%
    
    # ensure that this is arranged by date
    arrange(Date)
  
  # FOR NOW - assume the sex recorded in "fates" matches those in "covs"
  
  # attribute closest mass and HFL records to focal.id
  for (j in 1:nrow(focal.id)) {
    
    focal.id$Mass[j] <- focal.covs$Final.mass[which.min(abs(focal.id$start[j] - focal.covs$measure.date))] 
    
    focal.id$HFL[j] <- focal.covs$HFL[which.min(abs(focal.id$start[j] - focal.covs$measure.date))] 
    
  }
  
  # bind to master df
  fates.4 <- rbind(focal.id, fates.4)
  
}

# remove "start" and "end"
fates.5 <- fates.4 %>% dplyr::select(-c(start, end))

#_______________________________________________________________________________________________
# 7. Standardize and format covariates ----
#_______________________________________________________________________________________________

fates.5 <- fates.5 %>%
  
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
         Mass.1 = as.numeric(scale(Mass)),
         HFL.1 = as.numeric(scale(HFL))
         )

#_______________________________________________________________________________________________
# 8. Add treatment variable ----
#_______________________________________________________________________________________________

# here, we will use two indicator variables to estimate the effect of each treatment
# compared to control (i.e., the intercept)

# week 40 for 2A, 2B, 3A, and 3B
# week 41 for 1A, 1B, 4A, and 4B

fates.5$trt.ret <- NA
fates.5$trt.pil <- NA

# controls
fates.5$trt.ret[fates.5$Site %in% c("1C", "2C", "3C", "4C")] <- 0
fates.5$trt.pil[fates.5$Site %in% c("1C", "2C", "3C", "4C")] <- 0

# treatments in 2022 (akin to controls)
fates.5$trt.ret[fates.5$Site %in% c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B") &
                fates.5$year == 2022] <- 0
fates.5$trt.pil[fates.5$Site %in% c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B") &
                fates.5$year == 2022] <- 0

# 2 and 3 before week 40, 2023
fates.5$trt.ret[fates.5$Site %in% c("2A", "2B", "3A", "3B") &
                fates.5$week < 40 &
                fates.5$year == 2023] <- 0
fates.5$trt.pil[fates.5$Site %in% c("2A", "2B", "3A", "3B") &
                fates.5$week < 40 &
                fates.5$year == 2023] <- 0

# 1 and 4 before week 41, 2023
fates.5$trt.ret[fates.5$Site %in% c("1A", "1B", "4A", "4B") &
                fates.5$week < 41 &
                fates.5$year == 2023] <- 0
fates.5$trt.pil[fates.5$Site %in% c("1A", "1B", "4A", "4B") &
                fates.5$week < 41 &
                fates.5$year == 2023] <- 0

# rest of NAs for retention units
fates.5$trt.ret[is.na(fates.5$trt.ret) &
                      fates.5$Site %in% c("1A", "2B", "3B", "4A")] <- 1
fates.5$trt.pil[is.na(fates.5$trt.pil) &
                      fates.5$Site %in% c("1A", "2B", "3B", "4A")] <- 0

# rest of NAs for piling units
fates.5$trt.ret[is.na(fates.5$trt.ret) &
                      fates.5$Site %in% c("1B", "2A", "3A", "4B")] <- 0
fates.5$trt.pil[is.na(fates.5$trt.pil) &
                      fates.5$Site %in% c("1B", "2A", "3A", "4B")] <- 1

#_______________________________________________________________________________________________
# 9. Add cluster variable (will be an index) ----
#_______________________________________________________________________________________________

fates.5$cluster <- as.numeric(substr(fates.5$Site, 1, 1))

#_______________________________________________________________________________________________
# 10. Add principal component "body size" and body condition index ----
#_______________________________________________________________________________________________

# calculate first principal component axis
fates.5$PC1 <- (fates.5$Mass.1 + fates.5$HFL.1) / 2

# divide mass by HFL to create a body condition index uncorrelated with either variable
fates.5$BCI <- fates.5$Mass.1 / fates.5$HFL.1
fates.5$BCI.1 <- scale(fates.5$BCI)

#_______________________________________________________________________________________________
# 11. Write to csv ----
#_______________________________________________________________________________________________

write.csv(fates.5, "Cleaned data/fates_cleaned_08_30_2024.csv")
