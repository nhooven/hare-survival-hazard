# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04b - Model predictions (effect plots)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Dec 2023
# Date completed: 04 Mar 2024
# Date last modified: 04 Mar 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)
library(ggridges)        # ridgeline plot

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

load("RData - final/poisson_model.RData")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(m7))

#_______________________________________________________________________________________________
# 4. Visualize covariate effects and evaluate strength of evidence ----
#_______________________________________________________________________________________________
# 4a. All covariates (hazard ratios) - ridgeline plot ----
#_______________________________________________________________________________________________

# select parameters 
covariate.draws <- model.draws %>% dplyr::select(hr_sex:hr_bci)

# pivot
covariate.draws.long <- covariate.draws %>% 
  
  pivot_longer(cols = hr_sex:hr_bci) %>%
  
  # reorder factor
  mutate(name = factor(name,
                       levels = rev(c("hr_ret",
                                      "hr_pil",
                                      "hr_sex",
                                      "hr_mas",
                                      "hr_hfl",
                                      "hr_bci")),
                       labels = rev(c("Retention",
                                      "Piling",
                                      "Sex:M",
                                      "Body mass",
                                      "Hind foot length",
                                      "Body condition index"))))

# plot KDEs
ggplot(covariate.draws.long,
       aes(x = value,
           y = name)) +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(fill = "lightgray",
                      scale = 1,
                      rel_min_height = 0.00001) +
  
  # intercept at 1 (hazard ratio)
  geom_vline(xintercept = 1) +
  
  # theme
  theme(panel.grid = element_blank()) +
  
  # coordinates
  coord_cartesian(xlim = c(0.2, 2.5),
                  ylim = c(1.55, 6.25)) +
  
  # axis labels
  xlab("Hazard ratio") +
  ylab("Parameter")

#_______________________________________________________________________________________________
# 4b. All covariates (hazard ratios) - point estimates [medians] and CIs ----
#_______________________________________________________________________________________________

# here we'll use 95% and 50% CIs
sex.ci <- data.frame(med = median(covariate.draws$hr_sex),
                     lo.1 = quantile(covariate.draws$hr_sex, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_sex, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_sex, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_sex, prob = 0.75),
                     var = "Sex:M")

ret.ci <- data.frame(med = median(covariate.draws$hr_ret),
                     lo.1 = quantile(covariate.draws$hr_ret, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_ret, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_ret, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_ret, prob = 0.75),
                     var = "Retention")

pil.ci <- data.frame(med = median(covariate.draws$hr_pil),
                     lo.1 = quantile(covariate.draws$hr_pil, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_pil, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_pil, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_pil, prob = 0.75),
                     var = "Piling")

mas.ci <- data.frame(med = median(covariate.draws$hr_mas),
                     lo.1 = quantile(covariate.draws$hr_mas, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_mas, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_mas, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_mas, prob = 0.75),
                     var = "Body mass")

hfl.ci <- data.frame(med = median(covariate.draws$hr_hfl),
                     lo.1 = quantile(covariate.draws$hr_hfl, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_hfl, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_hfl, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_hfl, prob = 0.75),
                     var = "Hind foot length")

bci.ci <- data.frame(med = median(covariate.draws$hr_bci),
                     lo.1 = quantile(covariate.draws$hr_bci, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_bci, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_bci, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_bci, prob = 0.75),
                     var = "Body condition index")

# bind together
all.ci <- bind_rows(sex.ci, ret.ci, pil.ci, mas.ci, hfl.ci, bci.ci)

all.ci$var <- factor(all.ci$var,
                     levels = rev(c("Retention",
                                    "Piling",
                                    "Sex:M",
                                    "Body mass",
                                    "Hind foot length",
                                    "Body condition index")))

# plot
ggplot(data = all.ci,
       aes(x = med,
           y = var)) +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  
  # credible intervals
  # 95%
  geom_errorbarh(aes(xmin = lo.1,
                     xmax = up.1,
                     y = var),
                 height = 0,
                 linewidth = 2,
                 color = "lightgray") +
  # 50%
  geom_errorbarh(aes(xmin = lo.2,
                     xmax = up.2,
                     y = var),
                 height = 0,
                 linewidth = 2,
                 color = "darkgray") +
  
  # points
  geom_point(aes(x = med,
                 y = var),
             shape = 21,
             color = "black",
             fill = "white",
             size = 2.5,
             stroke = 1.1) +
  
  # axis titles
  xlab("Hazard ratio") +
  ylab("") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5))

#_______________________________________________________________________________________________
# 4c. Both density and point plots ----
#_______________________________________________________________________________________________

# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(data = covariate.draws.long,
                      aes(x = value,
                          y = name,
                          fill = name),
                      color = "#333333",
                      alpha = 0.5,
                      scale = 0.7,
                      rel_min_height = 0.005) +     # remove tails
  
  # credible intervals
  # 95%
  geom_errorbarh(data = all.ci,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = var,
                     color = var),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  # 50%
  geom_errorbarh(data = all.ci,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = var,
                     color = var),
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  
  # points
  geom_point(data = all.ci,
             aes(x = med,
                 y = var),
             shape = 21,
             color = "black",
             fill = "white",
             size = 2.5,
             stroke = 1.1,
             position = position_nudge(y = -0.1)) +
  
  # color
  scale_fill_manual(values = c("lightgray",
                               "lightgray",
                               "lightgray",
                               "lightgray",
                               "darkgreen",
                               "darkgreen")) +
  
  scale_color_manual(values = c("darkgray",
                                "darkgray",
                                "darkgray",
                                "darkgray",
                                "darkgreen",
                                "darkgreen")) +
  
  # vertical line at 1
  geom_vline(xintercept = 1) +
  
  # horizontal dashed line separating treatment and intrinsic variables
  geom_hline(yintercept = 4.55,
             linetype = "dashed",
             color = "darkgreen") +
  
  # axis titles
  xlab("Hazard ratio") +
  ylab("") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  
  # coordinates
  coord_cartesian(xlim = c(0.45, 2.65),
                  ylim = c(1.3, 6.3)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  # add text
  annotate(geom = "text",
           label = "Treatment",
           x = 2.4, 
           y = 6.65,
           color = "darkgreen",
           alpha = 0.5,
           size = 6) +
  
  annotate(geom = "text",
           label = "Intrinsic",
           x = 2.48, 
           y = 1,
           color = "darkgray",
           alpha = 0.75,
           size = 6)