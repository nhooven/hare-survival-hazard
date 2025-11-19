# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08 - Effect size plots
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2025 
# Date completed: 
# Date last modified: 19 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(ggridges)        # ridgeline plot
library(bayestestR)      # highest density posterior intervals

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- read.csv("Model outputs/model_1.csv")

#_______________________________________________________________________________________________
# 3. Visualize covariate effects and evaluate strength of evidence ----
#_______________________________________________________________________________________________
# 3a. Select hazard ratios and pivot_longer ----
#_______________________________________________________________________________________________

# select parameters 
covariate.draws <- model.fit.1 %>% dplyr::select(hr_bci,
                                                 hr_ret_total1,
                                                 hr_ret_total2,
                                                 hr_pil_total1,
                                                 hr_pil_total2)

# pivot
covariate.draws.long <- covariate.draws %>% 
  
  pivot_longer(cols = hr_bci:hr_pil_total2) %>%
  
  # reorder factor
  mutate(name = factor(name,
                       levels = rev(c("hr_pil_total1",
                                      "hr_pil_total2",
                                      "hr_ret_total1",
                                      "hr_ret_total2",
                                      "hr_bci")),
                       labels = rev(c("Piling (year 1)",
                                      "Piling (year 2)",
                                      "Retention (year 1)",
                                      "Retention (year 2)",
                                      "Body condition index"))))

#_______________________________________________________________________________________________
# 3b. Extract median and credible intervals ----
#_______________________________________________________________________________________________

# here we'll use 90% and 50% CIs
bci.ci <- data.frame(med = median(covariate.draws$hr_bci),
                     lo.1 = as.numeric(hdi(covariate.draws$hr_bci, ci = 0.90)[2]),
                     up.1 = as.numeric(hdi(covariate.draws$hr_bci, ci = 0.90)[3]),
                     lo.2 = as.numeric(hdi(covariate.draws$hr_bci, ci = 0.50)[2]),
                     up.2 = as.numeric(hdi(covariate.draws$hr_bci, ci = 0.50)[3]),
                     var = "Body condition index")

ret1.ci <- data.frame(med = median(covariate.draws$hr_ret_total1),
                      lo.1 = as.numeric(hdi(covariate.draws$hr_ret_total1, ci = 0.90)[2]),
                      up.1 = as.numeric(hdi(covariate.draws$hr_ret_total1, ci = 0.90)[3]),
                      lo.2 = as.numeric(hdi(covariate.draws$hr_ret_total1, ci = 0.50)[2]),
                      up.2 = as.numeric(hdi(covariate.draws$hr_ret_total1, ci = 0.50)[3]),
                      var = "Retention (year 1)")

ret2.ci <- data.frame(med = median(covariate.draws$hr_ret_total2),
                      lo.1 = as.numeric(hdi(covariate.draws$hr_ret_total2, ci = 0.90)[2]),
                      up.1 = as.numeric(hdi(covariate.draws$hr_ret_total2, ci = 0.90)[3]),
                      lo.2 = as.numeric(hdi(covariate.draws$hr_ret_total2, ci = 0.50)[2]),
                      up.2 = as.numeric(hdi(covariate.draws$hr_ret_total2, ci = 0.50)[3]),
                      var = "Retention (year 2)")

pil1.ci <- data.frame(med = median(covariate.draws$hr_pil_total1),
                      lo.1 = as.numeric(hdi(covariate.draws$hr_pil_total1, ci = 0.90)[2]),
                      up.1 = as.numeric(hdi(covariate.draws$hr_pil_total1, ci = 0.90)[3]),
                      lo.2 = as.numeric(hdi(covariate.draws$hr_pil_total1, ci = 0.50)[2]),
                      up.2 = as.numeric(hdi(covariate.draws$hr_pil_total1, ci = 0.50)[3]),
                      var = "Piling (year 1)")

pil2.ci <- data.frame(med = median(covariate.draws$hr_pil_total2),
                      lo.1 = as.numeric(hdi(covariate.draws$hr_pil_total2, ci = 0.90)[2]),
                      up.1 = as.numeric(hdi(covariate.draws$hr_pil_total2, ci = 0.90)[3]),
                      lo.2 = as.numeric(hdi(covariate.draws$hr_pil_total2, ci = 0.50)[2]),
                      up.2 = as.numeric(hdi(covariate.draws$hr_pil_total2, ci = 0.50)[3]),
                      var = "Piling (year 2)")



# bind together
all.ci <- bind_rows(bci.ci, ret1.ci, ret2.ci, pil1.ci, pil2.ci)

all.ci$var <- factor(all.ci$var,
                     levels = c("Body condition index",
                                "Retention (year 1)",
                                "Retention (year 2)",
                                "Piling (year 1)",
                                "Piling (year 2)"))

#_______________________________________________________________________________________________
# 3c. Both density and point plots ----
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
                      rel_min_height = 0.006) +     # remove tails
  
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
                               "darkgreen",
                               "darkgreen",
                               "darkgreen",
                               "darkgreen")) +
  
  scale_color_manual(values = c("darkgray",
                                "darkgreen",
                                "darkgreen",
                                "darkgreen",
                                "darkgreen")) +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "darkgray") +
  
  # axis titles
  xlab("Hazard ratio") +
  ylab("") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3)) +
  
  # coordinates
  coord_cartesian(xlim = c(0.35, 2.95),
                  ylim = c(1.3, 5)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black")) +
  
  # add text
  annotate(geom = "text",
           label = "Treatment",
           x = 2.3, 
           y = 5.35,
           color = "darkgreen",
           alpha = 0.5,
           size = 5) +
  
  annotate(geom = "text",
           label = "Intrinsic",
           x = 2.35, 
           y = 1.5,
           color = "darkgray",
           alpha = 0.75,
           size = 5)

