# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02b - Examine the censoring process
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 04 Apr 2024
# Date completed: 
# Date last modified: 04 Apr 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(lubridate)

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("Raw data/fates_04_04_2024.csv")

# initial data cleaning
fates.1 <- fates %>%
  
  # keep relevant columns
  dplyr::select(Site,
                Animal.ID,
                Sex,
                Ear.tag,
                Collar.type,
                Status,
                Capture.date,
                Estimated.event.date,
                Event.type,
                General.cause,
                Transmitter.lifetime) %>%
  
  # and observations
  filter(Event.type == "Censor" &
         General.cause %in% c("Unknown",
                              "Dead transmitter")) %>%
  
  # and convert to dates
  mutate(Capture.date = as.Date(mdy(Capture.date)),
         Estimated.event.date = as.Date(mdy(Estimated.event.date)))

#_______________________________________________________________________________________________
# 3. Visualize ----
#_______________________________________________________________________________________________

# density
ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ as.factor(Collar.type)) +
  
  geom_density(data = fates.1,
               aes(x = Transmitter.lifetime,
                   fill = as.factor(Collar.type)),
               color = "black",
               alpha = 0.5)  +
  
  theme(legend.position = "none",
        panel.grid = element_blank())

# histogram
ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ as.factor(Collar.type)) +
  
  geom_histogram(data = fates.1,
                 aes(x = Transmitter.lifetime,
                     fill = as.factor(Collar.type)),
                 color = "black",
                 alpha = 0.5,
                 binwidth = 28)  +
  
  theme(legend.position = "none",
        panel.grid = element_blank())
