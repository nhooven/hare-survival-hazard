# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 07 - Baseline hazard predictions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 
# Date last modified: 17 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- read.csv("Model outputs/model_1.csv")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

#_______________________________________________________________________________________________
# 3. Define spline ----
#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# weeks to predict on
weeks.pred <- seq(1, 52, length.out = 100)

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________________________
# 4. Prepare prediction data.frames and matrices ----
#_______________________________________________________________________________________________

# model fit as a df
model.fit.df <- as.data.frame(model.fit.1)

# subset model.draws with only the parameters we need
model.fit.spline <- model.fit.df %>%
  
  dplyr::select("a0.1.":"a0.4.",
                "w.1..1.":"w.4..9.")


# define 4 matrices with dimensions [n.draws, n.prediction points]
spline.pred.matrix.1 <- matrix(data = NA,
                               nrow = nrow(model.fit.spline),
                               ncol = nrow(basis.pred))

spline.pred.matrix.2 <- spline.pred.matrix.1
spline.pred.matrix.3 <- spline.pred.matrix.1
spline.pred.matrix.4 <- spline.pred.matrix.1

# here, we'll calculate pointwise predictions, one for each level (each sex/forest type)

#_______________________________________________________________________________________________
# 5. Define functions----

# This will be applied to each row of model.fit.spline, outputting to 
# spline.pred.array

# this function is exceedingly slow

#_______________________________________________________________________________________________
# 5a. Calculate splines first, then summarize ----
#_______________________________________________________________________________________________

spline_pred <- function(x,
                        sf) {
  
  # extract rowname for the array index [x, , ,]
  focal.index <- as.integer(rownames(x))
  
  # ensure we're working with the right s / ft
  focal.matrix <- case_when(sf == 1 ~ spline.pred.matrix.1,
                            sf == 2 ~ spline.pred.matrix.2,
                            sf == 3 ~ spline.pred.matrix.3,
                            sf == 4 ~ spline.pred.matrix.4)
  # predictions
  
  # convert coefficients to vector for multiplication
  w <- t(as.matrix(as.numeric(c(x[[paste0("w", ".", sf, "..1.")]],
                                x[[paste0("w", ".", sf, "..2.")]],
                                x[[paste0("w", ".", sf, "..3.")]],
                                x[[paste0("w", ".", sf, "..4.")]],
                                x[[paste0("w", ".", sf, "..5.")]],
                                x[[paste0("w", ".", sf, "..6.")]],
                                x[[paste0("w", ".", sf, "..7.")]],
                                x[[paste0("w", ".", sf, "..8.")]],
                                x[[paste0("w", ".", sf, "..9.")]]))))
  
  # multiply with 'sweep'
  w.by.b <- sweep(basis.pred, 2, w, `*`)
  
  # sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(w.by.b, 1, sum) 
  
  # add to the intercept (a0) and exponentiate to transform to hazard rate scale
  spline.pred <- exp(x[paste0("a0.", sf, ".")] + w.by.b.sum) 
  
  # and add into the array
  focal.matrix[focal.index, ] <- spline.pred
  
}

# apply function
spline.pred.1 <- apply(X = model.fit.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       sf = 1)

spline.pred.2 <- apply(X = model.fit.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       sf = 2)

spline.pred.3 <- apply(X = model.fit.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       sf = 3)

spline.pred.4 <- apply(X = model.fit.spline, 
                       MARGIN = 1, 
                       FUN = spline_pred,
                       sf = 4)

# calculate medians and upper/lower quantiles
spline.summary <- rbind(data.frame(x = weeks.pred,
                                   median = apply(spline.pred.1, 1, median),
                                   l90 = apply(spline.pred.1, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.1, 1, quantile, probs = 0.95),
                                   sf = 1),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.2, 1, median),
                                   l90 = apply(spline.pred.2, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.2, 1, quantile, probs = 0.95),
                                   sf = 2),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.3, 1, median),
                                   l90 = apply(spline.pred.3, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.3, 1, quantile, probs = 0.95),
                                   sf = 3),
                        data.frame(x = weeks.pred,
                                   median = apply(spline.pred.4, 1, median),
                                   l90 = apply(spline.pred.4, 1, quantile, probs = 0.05),
                                   u90 = apply(spline.pred.4, 1, quantile, probs = 0.95),
                                   sf = 4))

# sex and forest type variables
spline.summary <- spline.summary %>%
  
  mutate(sex = ifelse(sf %in% c(1, 2),
                      "F",
                      "M"),
         forest_type = ifelse(sf %in% c(1, 3),
                              "spruce-fir-lodgepole",
                              "xeric mixed conifer"))

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

#_______________________________________________________________________________________________
# 5a. Spline first / summarize later ----
#_______________________________________________________________________________________________

# plot
ggplot(data = spline.summary) +
  
  # facet
  facet_wrap(~ forest_type) +
  
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
                color = sex,
                linetype = sex),     
            linewidth = 1.25) +
  
  # ribbons for median spline prediction
  geom_ribbon(aes(x = x,
                  ymin = l90,
                  ymax = u90,
                  fill = sex),     
              alpha = 0.2) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 270,
                                   vjust = 0.25),
        legend.position = c(0.2, 0.7),
        legend.background = element_rect(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(fill = "white")) +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
                     labels = month.cutoffs$labels[c(2, 4, 6, 8, 10, 12)]) +
  
  scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09, 0.12)) +
  
  # axis labels
  ylab("Baseline hazard") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # colors
  scale_color_manual(values = c("pink2", "lightblue4")) +
  scale_fill_manual(values = c("pink2", "lightblue4"))

# 677 x 358

