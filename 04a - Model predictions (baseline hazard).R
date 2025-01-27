# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04a - Model predictions (baseline hazard)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Dec 2023
# Date completed: 04 Mar 2024
# Date last modified: 27 Jan 2024
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

load("RData - final/hazard_model.RData")

# fates data for week quantiles
fates <- read.csv("Cleaned data/fates_cleaned_01_02_2025.csv")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(hazard.model))

#_______________________________________________________________________________________________
# 4. Visualize baseline hazard function ----
#_______________________________________________________________________________________________
# 4a. Define basis functions ----
#_______________________________________________________________________________________________

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
weeks.pred <- seq(1, 52, length.out = 1000)

# define knots (quantile)
n.knots <- 7 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________________________
# 4b. Prepare prediction data.frames and matrices ----
#_______________________________________________________________________________________________

# subset model.draws with only the parameters we need
model.draws.spline <- model.draws %>%
  
  dplyr::select(a0_mean, 
                a0_sigma,
                a0_z.1, a0_z.2, a0_z.3, a0_z.4,
                w_pen.1.1:w_pen.4.7)

# define 4 matrices with dimensions [n.draws, n.prediction points]
spline.pred.matrix.1 <- matrix(data = NA,
                               nrow = nrow(model.draws.spline),
                               ncol = nrow(basis.pred))

spline.pred.matrix.2 <- spline.pred.matrix.1
spline.pred.matrix.3 <- spline.pred.matrix.1
spline.pred.matrix.4 <- spline.pred.matrix.1

# here, we'll calculate pointwise predictions, one for each level (each cluster)

#_______________________________________________________________________________________________
# 4c. Define a function ----

# This will be applied to each row of model.draws.spline, outputting to 
# spline.pred.array
#_______________________________________________________________________________________________

spline_pred <- function(x,
                        cluster,
                        ...) {
  
  # extract rowname for the array index [x, , ,]
  focal.index <- as.integer(rownames(x))
  
  # extract common parameters
  a0.mean <- x[["a0_mean"]]
  a0.sigma <- x[["a0_sigma"]]
  
  # ensure we're working with the right cluster
  focal.matrix <- case_when(cluster == 1 ~ spline.pred.matrix.1,
                            cluster == 2 ~ spline.pred.matrix.2,
                            cluster == 3 ~ spline.pred.matrix.3,
                            cluster == 4 ~ spline.pred.matrix.4)
  # predictions
  
  # convert coefficients to vector for multiplication
  w <- t(as.matrix(as.numeric(c(x[[paste0("w_pen.", cluster, ".1")]],
                                x[[paste0("w_pen.", cluster, ".2")]],
                                x[[paste0("w_pen.", cluster, ".3")]],
                                x[[paste0("w_pen.", cluster, ".4")]],
                                x[[paste0("w_pen.", cluster, ".5")]],
                                x[[paste0("w_pen.", cluster, ".6")]],
                                x[[paste0("w_pen.", cluster, ".7")]]))))
  
  # multiply with 'sweep'
  w.by.b <- sweep(basis.pred, 2, w, `*`)
  
  # sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(w.by.b, 1, sum) 
  
  # add to the intercept (a0) and exponentiate to transform to hazard rate scale
  spline.pred <- exp((a0.mean + a0.sigma * x[[paste0("a0_z.", cluster)]]) + w.by.b.sum) 
  
  # and add into the array
  focal.matrix[focal.index, ] <- spline.pred

}

# apply function
spline.pred.1 <- apply(X = model.draws.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       cluster = 1)

spline.pred.2 <- apply(X = model.draws.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       cluster = 2)

spline.pred.3 <- apply(X = model.draws.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       cluster = 3)

spline.pred.4 <- apply(X = model.draws.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       cluster = 4)

# calculate medians and upper/lower quantiles
spline.summary <- rbind(data.frame(x = weeks.pred,
                                   median = apply(spline.pred.1, 1, median),
                                   l90 = apply(spline.pred.1, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.1, 1, quantile, probs = 0.95),
                                   cluster = 1),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.2, 1, median),
                                   l90 = apply(spline.pred.2, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.2, 1, quantile, probs = 0.95),
                                   cluster = 2),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.3, 1, median),
                                   l90 = apply(spline.pred.3, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.3, 1, quantile, probs = 0.95),
                                   cluster = 3),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.4, 1, median),
                                   l90 = apply(spline.pred.4, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.4, 1, quantile, probs = 0.95),
                                   cluster = 4))

spline.summary$cluster <- factor(spline.summary$cluster,
                                 labels = c("Cluster 1",
                                            "Cluster 2",
                                            "Cluster 3",
                                            "Cluster 4"))
                           
#_______________________________________________________________________________________________
# 5. Plot median predictions and 90% quantiles ----
#_______________________________________________________________________________________________

# define month cutoffs
month.cutoffs <- data.frame(breaks = c(1, 5.57, 9.86, 
                                       14.29, 18.71, 22.71, 
                                       27.14, 31.43, 35.86, 
                                       40.14, 44.57, 48.86),
                            labels = c("Oct", "Nov", "Dec", 
                                       "Jan", "Feb", "Mar", 
                                       "Apr", "May", "Jun", 
                                       "Jul", "Aug", "Sep"))

# plot
ggplot(data = spline.summary) +
  
  # facet
  facet_wrap(~ cluster) +
  
  # white background
  theme_bw() +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "lightgray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "lightgray") +
  
  # line for median spline prediction
  geom_line(aes(x = x,
                y = median,
                color = cluster),     
            linewidth = 1.25) +
  
  # ribbons for median spline prediction
  geom_ribbon(aes(x = x,
                  ymin = l90,
                  ymax = u90,
                fill = cluster),     
            alpha = 0.2) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 270,
                                   vjust = 0),
        legend.position = "none",
        strip.text = element_text(hjust = 0)) +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
                     labels = month.cutoffs$labels[c(2, 4, 6, 8, 10, 12)]) +
  
  # axis labels
  ylab("Baseline hazard") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # colors
  scale_color_viridis_d(end = 0.9) +
  scale_fill_viridis_d(end = 0.9)

#_______________________________________________________________________________________________
# 6. Save image ----
#_______________________________________________________________________________________________

save.image(file = "RData - final/spline_pred.RData")
