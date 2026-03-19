# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 09c - Visualization - Landscape effects
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 09 Mar 2026 
# Date completed: 09 Mar 2026 
# Date last modified: 09 Mar 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# landscape metrics
lsm.demo <- read.csv("Cleaned data/lsm_demo.csv")

# fates data
fates <- read.csv("Cleaned data/fates_forModel.csv")

# model samples
model.fit <- readRDS("Model outputs/model_2.rds")

# convert to df
model.fit <- as.data.frame(do.call(rbind, model.fit))

# extract lsm effects only
model.fit <- model.fit %>% dplyr::select(all_of(c("b_dm", "b_open")))

#_______________________________________________________________________________________________
# 3. Data cleaning ----
#_______________________________________________________________________________________________
# 3a. LSM ----
#_______________________________________________________________________________________________

lsm.1 <- lsm.demo %>%
  
  # remove p.js
  dplyr::select(-p.js) %>%
  
  # new vars
  mutate(
    
    # standardize
    p.dm.s = (p.dm - mean(fates$p.dm)) / sd(fates$p.dm),
    p.o.s = (p.o - mean(fates$p.o)) / sd(fates$p.o),
    
    # standard site name
    site.s = c("1R", "1P", "1C", "2P", "2R", "2C", "3P", "3R", "3C", "4R", "4P", "4C")
    
  ) %>%
  
  # treatment
  mutate(
    
    trt = substr(site.s, 2, 2)
    
  )

#_______________________________________________________________________________________________
# 3b. Landscape effects ----

# here we want a 2D predictive surface, with posterior median and SD layers
# y axis will be DM, x-axis will be open

# let's start with 50 x 50

#_______________________________________________________________________________________________

ls_pred <- function (
    
  grid.size = 50
  
) {
  
  # generate (standardized) values to predict on
  vals.dm = seq(min(lsm.1$p.dm.s), max(lsm.1$p.dm.s), length.out = grid.size)
  vals.o = seq(min(lsm.1$p.o.s), max(lsm.1$p.o.s), length.out = grid.size)
  
  # expand grid (length = grid.size^2)
  vals.df <- expand.grid(vals.dm, vals.o)
  
  # loop through iterations
  hr.pred <- matrix(data = NA, 
                    nrow = nrow(model.fit),
                    ncol = grid.size * grid.size)
  
  for (i in 1:nrow(model.fit)) {
    
    focal.iter = as.matrix(model.fit[i, ])
    
    # calculate total effect
    hr.pred[i, ] = vals.df[ ,1] * focal.iter[1] + vals.df[ ,2] * focal.iter[2]
    
  }
  
  # exponentiate
  hr.pred.exp <- exp(hr.pred)
  
  # df for plotting
  hr.pred.df <- data.frame(
    
    p.dm.s = vals.df[ ,1],
    p.o.s = vals.df[ ,2],
    
    # calculate median and SD (by column)
    median = apply(hr.pred.exp, 2, median),
    sd = apply(hr.pred.exp, 2, sd)
    
  )
  
  return(hr.pred.df)
  
}

# apply function
ls.pred <- ls_pred(50)

#_______________________________________________________________________________________________
# 4. Plot ----

# breaks and labels
x.breaks = c(-2, -1, 0, 1)
y.breaks = c(-1, 0, 1, 2)

x.labels = round((x.breaks * sd(fates$p.o)) + mean(fates$p.o))
y.labels = round((y.breaks * sd(fates$p.dm)) + mean(fates$p.dm))

#_______________________________________________________________________________________________

# median
ggplot() +
  
  theme_classic() +
  
  # predictive surface
  geom_raster(
    
    data = ls.pred,
    aes(x = p.o.s,
        y = p.dm.s,
        fill = median)
    
  ) +
  
  # units
  geom_point(
    
    data = lsm.1,
    aes(x = p.o.s,
        y = p.dm.s)
    
  ) +

  scale_fill_gradient2(low = "purple",
                       high = "orange",
                       midpoint = 1) +
  
  # theme arguments
  theme(axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  
  # axis labels
  scale_x_continuous(breaks = x.breaks,
                     labels = x.labels) +
  
  scale_y_continuous(breaks = y.breaks,
                     labels = y.labels) +
  
  # labels
  xlab("% open") +
  ylab("% dense mature") +
  
  geom_text(data = lsm.1,
             aes(x = p.o.s,
                 y = p.dm.s,
                 label = site),
            nudge_x = 0.15)

# SD
ggplot() +
  
  theme_classic() +
  
  # predictive surface
  geom_raster(
    
    data = ls.pred,
    aes(x = p.o.s,
        y = p.dm.s,
        fill = sd)
    
  ) +
  
  # units
  geom_point(
    
    data = lsm.1,
    aes(x = p.o.s,
        y = p.dm.s),
    color = "white"
    
  ) +
  
  scale_fill_viridis_c() +
  
  # theme arguments
  theme(axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  
  # axis labels
  scale_x_continuous(breaks = x.breaks,
                     labels = x.labels) +
  
  scale_y_continuous(breaks = y.breaks,
                     labels = y.labels) +
  
  # labels
  xlab("% open") +
  ylab("% dense mature") +
  
  geom_text(data = lsm.1,
            aes(x = p.o.s,
                y = p.dm.s,
                label = site),
            nudge_x = 0.15)
