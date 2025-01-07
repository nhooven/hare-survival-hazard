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
# 2. Baseline hazard (old) ----
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

#_______________________________________________________________________________________________
# 3. Baseline hazard (spline with partial pooling) ----

# reference : https://www.tjmahr.com/random-effects-penalized-splines-same-thing/

# here we can learn a smoothing parameter, lambda, which can scale
# non-centered scaling factors for each w0 to control their smoothing level
# i.e., how shrunk they are toward zero.
# this constrains the w0 terms just like a skeptical prior BUT
# allows for partial pooling so they can learn from each other, rather
# than being fully independent

#_______________________________________________________________________________________________
# 3a. Toy example ----

library(mgcv)

#_______________________________________________________________________________________________

# create a synthetic dataset
weeks <- seq(1, 52, length.out = 100)

# define knots (quantile)
n.knots <- 15 + 1

knot.list <- quantile(weeks, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# creating the basis matrix (cyclic spline)
basis <- cSplineDes(weeks,
                    knots = knot.list,
                    ord = 4)                # cubic

# create some spline weights that are drawn from the same distribution,
# given some underlying factor (the non-centered scaling factor z)
z <- rnorm(n = n.knots - 1, mean = 0, sd = 1)

# define a sigma
sigma <- 0.25

# calculate spline weights
w0 <- sigma * z

# convert to vector for multiplication
w0.mat <- t(as.matrix(as.numeric(w0)))

# multiply with 'sweep'
w0.by.b <- sweep(basis, 2, w0, `*`)

# sum to create spline prediction
w0.by.b.sum <- apply(w0.by.b, 1, sum)

# create df to hold predictions
spline.pred.df <- data.frame(x = weeks,
                             y = w0.by.b.sum)

# plot
ggplot(spline.pred.df,
       aes(x = x,
           y = y)) +
  
  geom_point()

# what if we want to introduce a smoothing penalty?
# if lambda = 1, then there is no smoothing
# if lambda > 1, then smoothing is introduced
lambda <- 3

# calculate penalized weights
w0.penalized <- c(sigma * z[1] / lambda,
                  sigma * z[2] / lambda,
                  sigma * z[3] / lambda,
                  sigma * z[4] / lambda,
                  sigma * z[5] / lambda)

# convert to vector for multiplication
w0.penalized.mat <- t(as.matrix(as.numeric(w0.penalized)))

# multiply with 'sweep'
w0.penalized.by.b <- sweep(basis, 2, w0.penalized, `*`)

# sum to create spline prediction
w0.penalized.by.b.sum <- apply(w0.penalized.by.b, 1, sum)

spline.pred.df$y.penalized <- w0.penalized.by.b.sum

# plot
ggplot(spline.pred.df) +
  
  theme_bw() + 
  
  geom_point(aes(x = x,
                 y = y),
             color = "orange") +
  
  geom_point(aes(x = x,
                 y = y.penalized),
             color = "purple")
