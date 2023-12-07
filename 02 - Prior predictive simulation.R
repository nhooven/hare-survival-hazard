# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02 - Prior predictive simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 06 Dec 2023
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)

#_______________________________________________________________________________________________
# 2. Baseline hazard ----
#_______________________________________________________________________________________________

# Initial runs of the model with our data suggest that it struggles to estimate a constant
# baseline hazard - this makes sense since this should not be constant, but it will be 
# important to tune this based upon a sensible range of values and examine the MCMC
# consequences

# First: what would we expect weekly survival to be? Maybe we have more intuition about
# monthly survival.

# Let's say that for a given, average hare, monthly survival is 70%, with some small variability
surv.month <- rnorm(1000, 0.70, 0.1)

# now let's convert to weekly survival
surv.week <- surv.month ^ (1 / 4)

# and plot the distribution:
ggplot(data = as.data.frame(surv.week),
       aes(x = surv.week)) +
       geom_density(linewidth = 1.5) +
       theme_bw()

# now let's examine the consequences for the weekly baseline hazard rate:
haz.week <- -log(surv.week)

# and plot the distribution:
ggplot(data = as.data.frame(haz.week),
       aes(x = haz.week)) +
       geom_density(linewidth = 1.5) +
       theme_bw()

# this is pretty close to our initial posterior (mean = 0.17)
# now, let's parameterize a gamma distribution to approximate this
haz.week.gam <- rgamma(1000, 1, 1)

ggplot() +
       theme_bw() +
       geom_density(data = as.data.frame(haz.week),
                    aes(x = haz.week),
                    linewidth = 1.25,
                    color = "black") +
       geom_density(data = as.data.frame(haz.week.gam),
                    aes(x = haz.week.gam),
                    linewidth = 1.25,
                    color = "red")

# clearly the gamma(1, 1) prior is awful. We can do better and shrink this towards zero!
haz.week.gam <- rgamma(1000, 2, 20)

ggplot() +
  theme_bw() +
  geom_density(data = as.data.frame(haz.week),
               aes(x = haz.week),
               linewidth = 1.25,
               color = "black") +
  geom_density(data = as.data.frame(haz.week.gam),
               aes(x = haz.week.gam),
               linewidth = 1.25,
               color = "red")

# this looks better - let's give it a shot in the model!