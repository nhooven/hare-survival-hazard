# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04a - Model predictions (baseline hazard)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Dec 2023
# Date completed: 04 Mar 2024
# Date last modified: 28 Sep 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(mgcv)
library(ggridges)        # ridgeline plot

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

load("RData - final/poisson_model.RData")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(m3))

#_______________________________________________________________________________________________
# 4. Visualize time-varying baseline hazard function ----
#_______________________________________________________________________________________________
# 4a. Prepare predictions ----
#_______________________________________________________________________________________________

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
weeks.pred <- seq(1, 52, length.out = 100)           # 100 increments to predict on

basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

# we'll start by calculating predictions for every draw
spline.preds.df <- data.frame()

for (i in 1:nrow(model.draws)) {
  
  # subset spline parameters only
  spline.params <- model.draws[i , c(1, 7:11)]
  
  # convert to vector for multiplication
  w0 <- t(as.matrix(as.numeric(spline.params[ ,2:6])))
  
  # multiply with 'sweep'
  w0.by.b <- sweep(basis.pred, 2, w0, `*`)
  
  # sum to create "normal" scale spline prediction
  w0.by.b.sum <- apply(w0.by.b, 1, sum)
  
  # add to the intercept (a0) and exponentiate to transform to hazard rate scale
  spline.pred <- exp(spline.params$a0 + w0.by.b.sum)
  
  # create df to hold predictions
  spline.pred.df <- data.frame(x = weeks.pred,
                               y = spline.pred,
                               draw = i)
  
  # bind into master df
  spline.preds.df <- bind_rows(spline.preds.df,
                               spline.pred.df)
  
}

# summarize the predictions with quantiles
spline.preds.df.med <- spline.preds.df %>%
  
  group_by(x) %>%
  
  summarise(q = list(quantile(y, probs = c(0.025,
                                           0.05,
                                           0.10,
                                           0.125,
                                           0.50,
                                           0.875,
                                           0.90,
                                           0.95,
                                           0.975)))) %>% 
  
  unnest_wider(q)

#_______________________________________________________________________________________________
# 4b. Plot median predictions with all draws ----
#_______________________________________________________________________________________________

ggplot() +
  
  # white background
  theme_bw() +
  
  # line for median spline prediction
  geom_line(data = spline.preds.df.med,
            aes(x = x,
                y = `50%`),     
            linewidth = 1.25) +
  
  # lines for all predictions
  geom_line(data = spline.preds.df,
            aes(x = x,
                y = y,
                group = draw),
            alpha = 0.03) +
  
  # axis labels
  ylab("Baseline hazard") +
  xlab("Week of year")

#_______________________________________________________________________________________________
# 4c. Plot median predictions with credible intervals ----
#_______________________________________________________________________________________________

# pivot longer for plotting with a legend
spline.preds.df.med.long <- spline.preds.df.med %>%
  
  pivot_longer(cols = 2:ncol(spline.preds.df.med)) %>%
  
  mutate(level = ifelse(name %in% c("2.5%", "97.5%"),
                        "95%",
                        ifelse(name %in% c("5%", "95%"),
                               "90%",
                               ifelse(name %in% c("12.5%", "87.5%"),
                                      "75%",
                                      "50%"))))

# add repeating pattern for lower and upper intervals
spline.pred.ci <- spline.preds.df.med.long %>% 
  
  filter(level != "50%") %>%
  
  mutate(bound = rep_len(c("lo", "lo", "lo", "up", "up", "up"), n())) %>%
  
  dplyr::select(-name) %>%
  
  pivot_wider(names_from = c("bound")) %>%
  
  mutate(level = factor(level, levels = c("95%", "90%", "75%")))

# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # lines for all 4000 predictions
  geom_ribbon(data = spline.pred.ci,
              aes(x = x,
                  ymin = lo,
                  ymax = up,
                  alpha = level)) +
  
  # define alpha levels
  scale_alpha_discrete(range = c(0.05, 0.25)) +
  
  # line for median spline prediction
  geom_line(data = spline.preds.df.med,
            aes(x = x,
                y = `50%`),     
            linewidth = 1.25) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  # scale x axis by month
  scale_x_continuous(breaks = c(0, 4, 8, 12,
                                16, 20, 24, 28, 32, 
                                36, 40, 44, 48)) +
  
  # lines at the equinoxes and solstices
  geom_vline(xintercept = c(12, 25.5, 38.5, 51.5),
             alpha = 0.10,
             linewidth = 1.25,
             linetype = "dashed") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # axis labels
  ylab("Baseline hazard") +
  xlab("Week of year")
