# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 01 - Simulate time-to-event data and fit model
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 11 Apr 2023
# Date completed: 
# Date last modified: 25 Apr 2023
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)
library(survival)
library(rethinking)

#_______________________________________________________________________________________________
# 2. Simulate time-to-event data ----
#_______________________________________________________________________________________________

# we'll simulate data for the snow-on season (1 Nov-30 Apr)

#_______________________________________________________________________________________________
# 2a. Define covariates and simulated relationships ----
#_______________________________________________________________________________________________

# INTRINSIC
# SEX (F [0] or M [1]) - slightly higher hazard for males
b1 <- 0.25

# MASS - slightly lower hazard for heavier individuals
b2 <- -0.25

# SITE
# SITE - varying hazards per site
site.hazards <- data.frame("site" = c("C1", "T1", "C2", "T2", "C3", "T3"),
                           "site.hazard" = rnorm(6, 0, 0.05))

# TREATMENT (unthinned [0] or thinned [1]) - higher hazard for thinned
b3 <- 0.75

# TIME-VARYING
# SNOW (weekly snow depth) - slightly higher hazard for lower snow depths
b4 <- -0.25

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

# add baseline hazard rate 
# because S = exp(-blhaz), we need to simulate from a reasonable distribution
# let's say our 6-month survival follows a "tall" beta distribution
surv.rates <- rbeta(n = nrow(data.intrin), 8, 2)

data.intrin$bl.haz <- -log(surv.rates)

#_______________________________________________________________________________________________
# 2c. Define tracking period ----
#_______________________________________________________________________________________________

# ASSUME ALL ANIMALS WERE CAUGHT ON THE SAME DAY

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

# let's create 7 points (roughly translating to the first of each month)
snow.peaks <- data.frame(x = seq(1, 181, 30),
                         y = c(0, 40, 150, 125, 135, 90, 80))

ggplot(data = snow.peaks,
       aes(x = x,
           y = y)) +
  theme_bw() +
  geom_line(linetype = "dashed") +
  geom_point()

# create end points, slopes, and intercepts
for (i in 1:nrow(snow.peaks)) {
  
  snow.peaks$xend[i] <- snow.peaks$x[i + 1]
  
  snow.peaks$yend[i] <- snow.peaks$y[i + 1]
  
  snow.peaks$slope[i] <- (snow.peaks$yend[i] - snow.peaks$y[i])/
                         (snow.peaks$xend[i] - snow.peaks$x[i])
  
  snow.peaks$a[i] <- snow.peaks$y[i] - snow.peaks$slope[i]*snow.peaks$x[i]
  
}

# simulate data between points
sim.snow.pts <- list()

for (i in 1:6) {
  
  x.seq <- seq(snow.peaks$x[i], snow.peaks$xend[i], length = 100)
  
  sim.snow.pts[[i]] <- data.frame(x = x.seq,
                                  y = snow.peaks$a[i] + x.seq*snow.peaks$slope[i])
  
}

sim.snow.pts <- do.call(rbind, sim.snow.pts)

ggplot(data = snow.peaks,
       aes(x = x,
           y = y)) +
  theme_bw() +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_point(data = sim.snow.pts,
             aes(x = x,
                 y = y),
             color = "red")

# fit a GAM to these data
snow.model <- gam(y ~ s(x, k = 7), data = sim.snow.pts)

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
                y = predict),
            color = "red",
            linewidth = 1.25)

# calculate weekly averages
weekly.snow <- data.frame("day" = 1:181,
                          "week" = weeks,
                          "snow" = new.snow.peaks$predict) %>%
               group_by(week) %>%
               summarize(mean = mean(snow))

#_______________________________________________________________________________________________
# 2e. Expand grid to include all days ----
#_______________________________________________________________________________________________

all.indiv.days <- expand.grid("AnimalID" = data.intrin$AnimalID,
                              "day" = 1:181)

# add in week
days.weeks <- data.frame("day" = 1:181,
                         "week" = weeks)

all.indiv.days <- all.indiv.days %>% inner_join(days.weeks)

# add in time-varying covariate data
all.indiv.days <- all.indiv.days %>% inner_join(weekly.snow)

# join everything together
data.all <- data.intrin %>% inner_join(all.indiv.days) %>%
                            rename(snow = mean)

#_______________________________________________________________________________________________
# 2f. Simulate outcomes for each individual ----
#_______________________________________________________________________________________________

set.seed(6732)

final.data.all <- data.frame()

# loop through each individual
for (i in unique(data.all$AnimalID)) {
  
  # subset data per individual
  focal.id <- i
  
  focal.data <- data.all %>% filter(AnimalID == focal.id)
  
  # random censoring 
  # whether an animal is censored or not (Bernoulli trial, prob of success = 0.2)
  censor.yes <- rbinom(1, 1, 0.2)
  
  if (censor.yes == 1) {
    
    # if censored, draw day of censoring from a uniform distribution
    censor.day <- round(runif(1, 1, 181))
    
    focal.data <- focal.data %>% slice(1:(censor.day - 1)) %>%
                                 mutate(status = 0)
    
    # bind data into master df
    final.data.all <- bind_rows(final.data.all, focal.data)
    
  } else {
    
    # calculate hazard function for all individuals that can die
    hazard.data <- focal.data
    
    # determine index values
    index.sex <- ifelse(hazard.data$Sex[1] == "F", 0, 1)
    
    index.thin <- ifelse(hazard.data$Treatment[1] == "unthinned", 0, 1)
    
    site.effect <- site.hazards$site.hazard[site.hazards$site == hazard.data$Site[1]]
    
    # calculate hazard
    hazard.data <- hazard.data %>% mutate(hazard = bl.haz *              # baseline hazard
                                                   exp(
                                                       site.effect +
                                                       b1 * index.sex +
                                                       b2 * Mass +
                                                       b3 * index.thin +
                                                       b4 * snow
                                                   )) %>%
                                   mutate(log.hazard = log(hazard))

    # death as a Bernoulli trial
    hazard.data$status <- rbinom(nrow(hazard.data), 1, prob = plogis(hazard.data$log.hazard))
    
    # keep only rows before the 1
    first.fate <- which(hazard.data$status == 1, arr.ind = TRUE)[1]
    
    if (is.na(first.fate) == FALSE) {
      
      hazard.data <- hazard.data %>% slice(1:first.fate)
      
    }
    
    # bind data into master df
    final.data.all <- bind_rows(final.data.all, hazard.data)
    
  }
  
}

# keep only columns we need for modeling
final.data.all.1 <- final.data.all %>% dplyr::select(AnimalID, 
                                                     Sex, 
                                                     Mass, 
                                                     Site, 
                                                     Treatment, 
                                                     snow, 
                                                     day, 
                                                     week, 
                                                     status)

# aggregate to weekly survival data
final.data.all.2 <- data.frame()
  
for (i in unique(final.data.all.1$AnimalID)) {
  
  # subset data per individual
  focal.id <- i
  
  focal.data <- final.data.all.1 %>% filter(AnimalID == focal.id)
  
  for (x in unique(focal.data$week)) {
    
    # subset data by week
    focal.week <- x
    
    focal.week.data <- focal.data %>% filter(week == focal.week)
    
    # new focal week data
    focal.week.data.1 <- focal.week.data %>% slice(1) %>%
                                             dplyr::select(AnimalID, 
                                                           Sex, 
                                                           Mass, 
                                                           Site, 
                                                           Treatment, 
                                                           snow, 
                                                           week)
    
    focal.week.data.1$status <- max(focal.week.data$status)
    
    # bind into master df
    final.data.all.2 <- bind_rows(final.data.all.2, focal.week.data.1)
    
    }
  
}

# add start and end week columns
final.data.all.3 <- final.data.all.2 %>% mutate(start.week = as.integer(week) - 1,
                                                end.week = week) %>%
                                         dplyr::select(-week)

#_______________________________________________________________________________________________
# 3. Summarize and visualize data ----
#_______________________________________________________________________________________________

# median n of weeks survived
median(as.integer(final.data.all.3$end.week))

# mean n of weeks survived
min(as.integer(final.data.all.3$end.week))

# max of weeks survived
max(as.integer(final.data.all.3$end.week))

# distribution of weeks survived
ggplot(data = final.data.all.3,
       aes(x = as.integer(end.week))) +
  theme_bw() +
  geom_density(linewidth = 1.25) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(1, 26, 2)) +
  xlab("Week of death")

# distribution by Sex
ggplot(data = final.data.all.3,
       aes(x = as.integer(end.week),
           color = Sex,
           fill = Sex)) +
  theme_bw() +
  geom_density(linewidth = 1.25,
               alpha = 0.25) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.8)) +
  scale_x_continuous(breaks = seq(1, 26, 2)) +
  xlab("Week of death")

# distribution by Site
ggplot(data = final.data.all.3,
       aes(x = as.integer(end.week),
           color = Site,
           fill = Site)) +
  theme_bw() +
  geom_density(linewidth = 1.25,
               alpha = 0.25) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.8)) +
  scale_x_continuous(breaks = seq(1, 26, 2)) +
  xlab("Week of death")

#_______________________________________________________________________________________________
# 4. Fit Cox regression in Stan ----
#_______________________________________________________________________________________________

# scale continuous variables
final.data.all.3$Mass.s <- scale(final.data.all.3$Mass)
final.data.all.3$snow.s <- scale(final.data.all.3$snow)

# estimate baseline hazard
final.data.baseline <- list(status = final.data.all.3$status,
                            start_week = final.data.all.3$start.week)

m1 <- ulam(
  
  alist(
    
    status ~ dpois(mu),        
    mu <- bl,
    
    # priors
    bl ~ dexp(1)
    
  ), data = final.data.baseline, 
     chains = 4, 
     cores = 4
)

precis(m1)

# estimate covariate effects - mass
final.data.cov <- list(status = final.data.all.3$status,
                       start_week = final.data.all.3$start.week,
                       mass = final.data.all.3$Mass.s)

m2 <- ulam(
  
  alist(
    
    status ~ dpois(mu),        
    mu <- bl * exp(b1*mass),
    
    # priors
    bl ~ dexp(1),
    b1 ~ dnorm(0, 1)
    
  ), data = final.data.cov, 
  chains = 4, 
  cores = 4
)

precis(m2)
traceplot(m2)
