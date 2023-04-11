# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Simulate time-to-event data and fit model
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 11 Apr 2023
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Simulate time-to-event data ----
#_______________________________________________________________________________________________

# we'll simulate data for the snow-on season (1 Nov-30 Apr)

#_______________________________________________________________________________________________
# 2a. Define covariates and simulated relationships ----
#_______________________________________________________________________________________________

# INTRINSIC
# SEX (F or M) - slightly higher hazard for males
# MASS - slightly lower hazard for heavier individuals

# SITE
# SITE - varying hazards per site
# TREATMENT (unthinned or thinned) - higher hazard for thinned

# TIME-VARYING
# SNOW (weekly snow depth) - slightly higher hazard for lower snow depths
# TEMP (weekly air temperature) - slightly higher hazard for lower temperatures

#_______________________________________________________________________________________________
# 2b. Simulate intrinsic and site covariates for 100 individuals ----
#_______________________________________________________________________________________________

n.indiv <- 100

# add sex variable
data.intrin <- tibble("AnimalID" = as.factor(1:n.indiv),      
                      "Sex" = rep(c("F", "M"), 50),      # equal sex ratio
                      "Mass" = NA)       

# make females slightly heavier than males to mirror reality
for (i in 1:nrow(data.intrin)) {
  
  data.intrin$Mass[i] <- ifelse("Sex" == "F",
                                      rnorm(n = length(n.indiv),
                                            mean = 1.2,
                                            sd = 0.25),
                                      rnorm(n = length(n.indiv),
                                            mean = 1.0,
                                            sd = 0.25))
  
}

# check that these are reasonable
range(data.intrin$Mass)

# add in study sites (one of 6)
sites <- c("C1", "T1", "C2", "T2", "C3", "T3")

data.intrin$Site <- as.factor(sample(sites, size = 100, replace = TRUE))

# add treatment
data.intrin$Treatment <- as.factor(substr(data.intrin$Site, 1, 1))

#_______________________________________________________________________________________________
# 2c. Define tracking period ----
#_______________________________________________________________________________________________

# how many days?
total.days <- 30 + 31 + 31 + 28 + 31 + 30

# how many weeks?
weeks <- cut(1:total.days, 
             breaks = seq(0, 181 + 7, 7),
             labels = c(1:26))

total.weeks <- length(unique(weeks))

#_______________________________________________________________________________________________
# 2d. Simulate daily weather data ----
#_______________________________________________________________________________________________

# https://www.r-bloggers.com/2022/08/simulating-data-from-a-non-linear-function-by-specifying-a-handful-of-points/

# snow depth (in cm) is a nonlinear function of day, with a peak on Jan 1 (day 62)

# let's create 6 points (roughly translating to the first of each month)
snow.peaks <- data.frame(x = seq(1, 151, 30),
                         y = c(0, 40, 150, 125, 135, 90))

ggplot(data = snow.peaks,
       aes(x = x,
           y = y)) +
  theme_bw() +
  geom_line(linetype = "dashed") +
  geom_point()

# fit a GAM to these data
snow.model <- gam(y ~ s(x, k = 7), data = snow.peaks)

# and predict on a new data frame
new.snow.peaks <- data.frame(x = 1:181)

new.snow.peaks$predict <- predict(snow.model, new.snow.peaks)

ggplot(data = snow.peaks,
       aes(x = x,
           y = y)) +
  theme_bw() +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_line(data = new.snow.peaks,
            aes(x = x,
                y = predict))

#_______________________________________________________________________________________________
# 2e. Expand grid to include all days ----
#_______________________________________________________________________________________________

