# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Clean fates data
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2023
# Date completed: 27 Nov 2023
# Date last modified: 01 Jan 2024 
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

fates <- read.csv("fates_covariates_11_20_2023.csv")

#_______________________________________________________________________________________________
# 3. Select only columns/rows we need ----
#_______________________________________________________________________________________________

fates.1 <- fates %>% 
  
  # remove hare with collar 150.533 for now
  filter(Notes != "Unknown hare") %>%
  
  # select relevant variables
  dplyr::select(
  
    Site,
    Animal.ID,
    Sex,
    Ear.tag,
    Collar.type,
    Mass,
    HFL,
    Status,
    Capture.date,
    Estimated.event.date,
    Event.type,
    Include
    
  ) %>%
  
  # filter out those not to include (e.g., capture-related)
  filter(Include == "Y") %>%
  
  # remove "Include" variables
  dplyr::select(-Include) 
  
#_______________________________________________________________________________________________
# 4. Impute missing covariate data ----
#_______________________________________________________________________________________________

# determine mean mass (rounded)
mean.mass <- round(mean(fates.1$Mass, na.rm = T), digits = 2)
mean.mass.F <- round(mean(fates.1$Mass[fates.1$Sex == "F"], na.rm = T), digits = 2)
mean.mass.M <- round(mean(fates.1$Mass[fates.1$Sex == "M"], na.rm = T), digits = 2)

# determine mean HFL (rounded)
mean.hfl <- round(mean(fates.1$HFL, na.rm = T), digits = 1)
mean.hfl.F <- round(mean(fates.1$HFL[fates.1$Sex == "F"], na.rm = T), digits = 1)
mean.hfl.M <- round(mean(fates.1$HFL[fates.1$Sex == "M"], na.rm = T), digits = 1)

# add mean values to each missing value
# mass
fates.1$Mass[fates.1$Ear.tag == 399] <- mean.mass.F
fates.1$Mass[fates.1$Ear.tag == 986] <- mean.mass.F
fates.1$Mass[fates.1$Ear.tag == 989] <- mean.mass.M
fates.1$Mass[fates.1$Ear.tag == 1883] <- mean.mass.M

# HFL
fates.1$HFL[fates.1$Ear.tag == 393] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 399] <- mean.hfl.F
fates.1$HFL[fates.1$Ear.tag == 456] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 478] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 483] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 527] <- mean.hfl.F
fates.1$HFL[fates.1$Ear.tag == 986] <- mean.hfl.F
fates.1$HFL[fates.1$Ear.tag == 989] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 1061] <- mean.hfl.F
fates.1$HFL[fates.1$Ear.tag == 1062] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 1730] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 1737] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == 1883] <- mean.hfl.M
fates.1$HFL[fates.1$Ear.tag == "399 / 1063"] <- mean.hfl.F

#_______________________________________________________________________________________________
# 5. Format dates correctly ----
#_______________________________________________________________________________________________

fates.1$Capture.date <- as.Date(mdy(fates.1$Capture.date))
fates.1$Estimated.event.date <- as.Date(mdy(fates.1$Estimated.event.date))

# replace NAs as most recently updated date
fates.1$Estimated.event.date[is.na(fates.1$Estimated.event.date)] <- as.Date("2023-11-20")

#_______________________________________________________________________________________________
# 6. Split into time periods by week ----
#_______________________________________________________________________________________________

# here, we will use the survSplit function in the "survival" package to split the dataset by week

# create numeric start and end variables
# define origin for our analysis (2022-01-01)
date.origin <- as.numeric(as.Date("2022-01-01"))

fates.1$start <- as.numeric(fates.1$Capture.date) - date.origin
fates.1$end <- as.numeric(fates.1$Estimated.event.date) - date.origin

# create vector of time points
cut.points <- seq(1,
                  as.numeric(as.Date("2023-11-20")) - date.origin,
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

# keep only columns we need
fates.3 <- fates.2 %>%
  
  # select
  dplyr::select(Site,
                Animal.ID,
                Ear.tag,
                Collar.type,
                Sex,
                Mass,
                HFL,
                Event.type,
                start,
                end,
                status.num)

# create numeric week variable (back-transform)
fates.3$week <- as.numeric(format(as.Date(fates.3$end + date.origin, 
                                          origin = "1970-01-01"), 
                                  "%W"))

# change 0 weeks to 52
fates.3$week[fates.3$week == 0] <- 52

# add year variable
fates.3$year <- ifelse(fates.3$end < 366,
                       2022,
                       2023)

#_______________________________________________________________________________________________
# 7. Standardize covariates ----
#_______________________________________________________________________________________________

# use an indicator variable for sex (0 = F [intercept])
fates.3$Sex.1 <- ifelse(fates.3$Sex == "F",
                        0,
                        1)

fates.3$Mass.1 <- as.numeric(scale(fates.3$Mass))

fates.3$HFL.1 <- as.numeric(scale(fates.3$HFL))

#_______________________________________________________________________________________________
# 8. Add treatment variable ----
#_______________________________________________________________________________________________

# here, we will use two indicator variables to estimate the effect of each treatment
# compared to control (i.e., the intercept)

# week 40 for 2A, 2B, 3A, and 3B
# week 41 for 1A, 1B, 4A, and 4B

fates.3$Treatment.Retention <- NA
fates.3$Treatment.Piling <- NA

# controls
fates.3$Treatment.Retention[fates.3$Site %in% c("1C", "2C", "3C", "4C")] <- 0
fates.3$Treatment.Piling[fates.3$Site %in% c("1C", "2C", "3C", "4C")] <- 0

# treatments in 2022 (akin to controls)
fates.3$Treatment.Retention[fates.3$Site %in% c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B") &
                            fates.3$year == 2022] <- 0
fates.3$Treatment.Piling[fates.3$Site %in% c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B") &
                         fates.3$year == 2022] <- 0

# 2 and 3 before week 40, 2023
fates.3$Treatment.Retention[fates.3$Site %in% c("2A", "2B", "3A", "3B") &
                            fates.3$week < 40 &
                            fates.3$year == 2023] <- 0
fates.3$Treatment.Piling[fates.3$Site %in% c("2A", "2B", "3A", "3B") &
                         fates.3$week < 40 &
                         fates.3$year == 2023] <- 0

# 1 and 4 before week 41, 2023
fates.3$Treatment.Retention[fates.3$Site %in% c("1A", "1B", "4A", "4B") &
                            fates.3$week < 41 &
                            fates.3$year == 2023] <- 0
fates.3$Treatment.Piling[fates.3$Site %in% c("1A", "1B", "4A", "4B") &
                         fates.3$week < 41 &
                         fates.3$year == 2023] <- 0

# rest of NAs for retention units
fates.3$Treatment.Retention[is.na(fates.3$Treatment.Retention) &
                            fates.3$Site %in% c("1A", "2B", "3B", "4A")] <- 1
fates.3$Treatment.Piling[is.na(fates.3$Treatment.Piling) &
                         fates.3$Site %in% c("1A", "2B", "3B", "4A")] <- 0

# rest of NAs for piling units
fates.3$Treatment.Retention[is.na(fates.3$Treatment.Retention) &
                            fates.3$Site %in% c("1B", "2A", "3A", "4B")] <- 0
fates.3$Treatment.Piling[is.na(fates.3$Treatment.Piling) &
                         fates.3$Site %in% c("1B", "2A", "3A", "4B")] <- 1

#_______________________________________________________________________________________________
# 9. Add cluster variable (will be an index) ----
#_______________________________________________________________________________________________

fates.3$cluster <- as.numeric(substr(fates.3$Site, 1, 1))

#_______________________________________________________________________________________________
# 10. Write to csv ----
#_______________________________________________________________________________________________

write.csv(fates.3, "fates_cleaned.csv")
                         
