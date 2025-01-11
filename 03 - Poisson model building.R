# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03a - Poisson model building
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 29 Dec 2023
# Date last modified: 10 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)         # construct basic functions
library(mgcv)            # cyclic splines

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_cleaned_01_02_2025.csv")

# keep only columns we need for modeling
# spline hazard only
fates.1 <- fates %>% dplyr::select(cluster,
                                   Site,
                                   Ear.tag,
                                   week,
                                   year,
                                   mort.pred,
                                   mort.oth,
                                   cens)

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 5 + 1

knot.list <- quantile(fates.1$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# set first knot to zero
knot.list[1] <- 0

# creating the basis matrix (cyclic spline)
basis <- cSplineDes(fates.1$week,
                    knots = knot.list,
                    ord = 4)                # cubic

# plot basis functions
basis.plot <- tibble(week = fates.1$week,
                     b1 = basis[ , 1],
                     b2 = basis[ , 2],
                     b3 = basis[ , 3],
                     b4 = basis[ , 4],
                     b5 = basis[ , 5]) %>%
  pivot_longer(cols = b1:b5)

# plot
ggplot(data = fates.1,
       aes(x = week)) +
  
  theme_bw() +
  
  geom_line(data = basis.plot,
            aes(x = week,
                y = value,
                group = name),
            linewidth = 1.5,
            alpha = 0.25)

#_______________________________________________________________________________________________
# 4. Build data list ----
#_______________________________________________________________________________________________

# splines
fates.stan.1 <- list(N = nrow(fates.1),
                     y_mort_pred = fates.1$mort.pred,
                     y_mort_oth = fates.1$mort.oth,
                     y_cens = fates.1$cens,
                     t = fates.1$week,
                     clust = fates.1$cluster,
                     nclust = 4,
                     basis = basis,
                     nbasis = ncol(basis))

#_______________________________________________________________________________________________
# 5. Testing and tuning the spline ----
#_______________________________________________________________________________________________
# 5a. Plot predicted spline function ----
#_______________________________________________________________________________________________

plot_spline <- function (weeks.pred,
                         spline.intercept,
                         spline.coefs) {
  
  # prediction basis function
  basis.pred <- cSplineDes(weeks.pred,
                           knots = knot.list,
                           ord = 4)
  
  # convert to vector for multiplication
  w0.penalized.mat <- t(as.matrix(as.numeric(spline.coefs)))
  
  # multiply with 'sweep'
  w0.penalized.by.b <- sweep(basis.pred, 2, w0.penalized.mat, `*`)
  
  # sum to create spline prediction
  w0.penalized.by.b.sum <- rowSums(w0.penalized.by.b)
  
  # add to intercept and exponentiate
  spline.pred <- exp(spline.intercept + w0.penalized.by.b.sum)
  
  # data.frame
  spline.pred.df <- data.frame(x = weeks.pred,
                               y = spline.pred)
  
  # plot
  ggplot(spline.pred.df,
         aes(x = x,
             y = y)) +
    
    theme_bw() + 
    
    geom_point(color = "purple") 
  
}

#_______________________________________________________________________________________________
# 5b. Run model ----
#_______________________________________________________________________________________________

m.spline <- rstan::stan(
  
  file = "model_spline_test.stan",
  data = fates.stan.1,
  chains = 1,
  warmup = 1000,
  iter = 2000
  
)

print(m.spline)

plot(m.spline, pars = "zw_pred")

# plot spline
plot_spline(seq(1, 52, length.out = 100),
            -3.86,
            c(0.20, -0.65, -0.16, -0.68, -0.31))



