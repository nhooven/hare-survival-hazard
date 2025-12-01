# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08 - Effect size plots
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2025 
# Date completed: 01 Dec 2025
# Date last modified: 01 Dec 2025
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
model.fit.2 <- read.csv("Model outputs/model_2.csv")
model.fit.3 <- read.csv("Model outputs/model_3.csv")

#_______________________________________________________________________________________________
# 3. Visualize covariate effects and evaluate strength of evidence ----
#_______________________________________________________________________________________________
# 3a. Select hazard ratios and pivot_longer ----
#_______________________________________________________________________________________________

# select parameters 
covariate.draws.1 <- model.fit.1 %>% dplyr::select(hr_bci,
                                                   hr_bci_study_week,
                                                   hr_ret_total1,
                                                   hr_ret_total2,
                                                   hr_pil_total1,
                                                   hr_pil_total2)

covariate.draws.2 <- model.fit.2 %>% dplyr::select(hr_bci,
                                                   hr_bci_study_week,
                                                   hr_ret_total1,
                                                   hr_ret_total2,
                                                   hr_pil_total1,
                                                   hr_pil_total2)

covariate.draws.3 <- model.fit.3 %>% dplyr::select(hr_bci,
                                                   hr_bci_study_week,
                                                   hr_ret_total1,
                                                   hr_ret_total2,
                                                   hr_pil_total1,
                                                   hr_pil_total2)

# pivot
covariate.draws.long.1 <- covariate.draws.1 %>% 
  
  pivot_longer(cols = hr_bci:hr_pil_total2) %>%
  
  # reorder factor
  mutate(name = factor(name,
                       levels = rev(c("hr_pil_total1",
                                      "hr_pil_total2",
                                      "hr_ret_total1",
                                      "hr_ret_total2",
                                      "hr_bci",
                                      "hr_bci_study_week")),
                       labels = rev(c("Piling (year 1)",
                                      "Piling (year 2)",
                                      "Retention (year 1)",
                                      "Retention (year 2)",
                                      "Body condition index",
                                      "Body condition index * t"))),
         scen = "1")

# model 2
covariate.draws.long.2 <- covariate.draws.2 %>% 
  
  pivot_longer(cols = hr_bci:hr_pil_total2) %>%
  
  # reorder factor
  mutate(name = factor(name,
                       levels = rev(c("hr_pil_total1",
                                      "hr_pil_total2",
                                      "hr_ret_total1",
                                      "hr_ret_total2",
                                      "hr_bci",
                                      "hr_bci_study_week")),
                       labels = rev(c("Piling (year 1)",
                                      "Piling (year 2)",
                                      "Retention (year 1)",
                                      "Retention (year 2)",
                                      "Body condition index",
                                      "Body condition index * t"))),
         scen = "2")

# model 3
covariate.draws.long.3 <- covariate.draws.3 %>% 
  
  pivot_longer(cols = hr_bci:hr_pil_total2) %>%
  
  # reorder factor
  mutate(name = factor(name,
                       levels = rev(c("hr_pil_total1",
                                      "hr_pil_total2",
                                      "hr_ret_total1",
                                      "hr_ret_total2",
                                      "hr_bci",
                                      "hr_bci_study_week")),
                       labels = rev(c("Piling (year 1)",
                                      "Piling (year 2)",
                                      "Retention (year 1)",
                                      "Retention (year 2)",
                                      "Body condition index",
                                      "Body condition index * t"))),
         scen = "3")

#_______________________________________________________________________________________________
# 3b. Extract median and credible intervals ----

# function
extract_med_ci <- function (x) {
  
  # here we'll use 90% and 50% CIs
  bci.ci <- data.frame(med = median(x$hr_bci),
                       lo.1 = as.numeric(hdi(x$hr_bci, ci = 0.90)[2]),
                       up.1 = as.numeric(hdi(x$hr_bci, ci = 0.90)[3]),
                       lo.2 = as.numeric(hdi(x$hr_bci, ci = 0.50)[2]),
                       up.2 = as.numeric(hdi(x$hr_bci, ci = 0.50)[3]),
                       var = "Body condition index")
  
  bci.t.ci <- data.frame(med = median(x$hr_bci_study_week),
                         lo.1 = as.numeric(hdi(x$hr_bci_study_week, ci = 0.90)[2]),
                         up.1 = as.numeric(hdi(x$hr_bci_study_week, ci = 0.90)[3]),
                         lo.2 = as.numeric(hdi(x$hr_bci_study_week, ci = 0.50)[2]),
                         up.2 = as.numeric(hdi(x$hr_bci_study_week, ci = 0.50)[3]),
                         var = "Body condition index * t")
  
  ret1.ci <- data.frame(med = median(x$hr_ret_total1),
                        lo.1 = as.numeric(hdi(x$hr_ret_total1, ci = 0.90)[2]),
                        up.1 = as.numeric(hdi(x$hr_ret_total1, ci = 0.90)[3]),
                        lo.2 = as.numeric(hdi(x$hr_ret_total1, ci = 0.50)[2]),
                        up.2 = as.numeric(hdi(x$hr_ret_total1, ci = 0.50)[3]),
                        var = "Retention (year 1)")
  
  ret2.ci <- data.frame(med = median(x$hr_ret_total2),
                        lo.1 = as.numeric(hdi(x$hr_ret_total2, ci = 0.90)[2]),
                        up.1 = as.numeric(hdi(x$hr_ret_total2, ci = 0.90)[3]),
                        lo.2 = as.numeric(hdi(x$hr_ret_total2, ci = 0.50)[2]),
                        up.2 = as.numeric(hdi(x$hr_ret_total2, ci = 0.50)[3]),
                        var = "Retention (year 2)")
  
  pil1.ci <- data.frame(med = median(x$hr_pil_total1),
                        lo.1 = as.numeric(hdi(x$hr_pil_total1, ci = 0.90)[2]),
                        up.1 = as.numeric(hdi(x$hr_pil_total1, ci = 0.90)[3]),
                        lo.2 = as.numeric(hdi(x$hr_pil_total1, ci = 0.50)[2]),
                        up.2 = as.numeric(hdi(x$hr_pil_total1, ci = 0.50)[3]),
                        var = "Piling (year 1)")
  
  pil2.ci <- data.frame(med = median(x$hr_pil_total2),
                        lo.1 = as.numeric(hdi(x$hr_pil_total2, ci = 0.90)[2]),
                        up.1 = as.numeric(hdi(x$hr_pil_total2, ci = 0.90)[3]),
                        lo.2 = as.numeric(hdi(x$hr_pil_total2, ci = 0.50)[2]),
                        up.2 = as.numeric(hdi(x$hr_pil_total2, ci = 0.50)[3]),
                        var = "Piling (year 2)")
  
  
  
  # bind together
  all.ci <- bind_rows(bci.ci, bci.t.ci, ret1.ci, ret2.ci, pil1.ci, pil2.ci)
  
  all.ci$var <- factor(all.ci$var,
                       levels = c("Body condition index",
                                  "Body condition index * t", 
                                  "Retention (year 1)",
                                  "Retention (year 2)",
                                  "Piling (year 1)",
                                  "Piling (year 2)"))
  
  return (all.ci)
  
}

#_______________________________________________________________________________________________

# use function
all.ci.1 <- extract_med_ci(covariate.draws.1)
all.ci.2 <- extract_med_ci(covariate.draws.2)
all.ci.3 <- extract_med_ci(covariate.draws.3)

# add scen
all.ci.1$scen = "1"
all.ci.2$scen = "2"
all.ci.3$scen = "3"

#_______________________________________________________________________________________________
# 4. By scenario ----
#_______________________________________________________________________________________________

# model 1
# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(data = covariate.draws.long.1,
                      aes(x = value,
                          y = name,
                          fill = name),
                      color = "#333333",
                      alpha = 0.5,
                      scale = 0.7,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  # 95%
  geom_errorbarh(data = all.ci.1,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = var,
                     color = var),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  # 50%
  geom_errorbarh(data = all.ci.1,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = var,
                     color = var),
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  
  # points
  geom_point(data = all.ci.1,
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
                               "darkgreen",
                               "darkgreen",
                               "darkgreen",
                               "darkgreen")) +
  
  scale_color_manual(values = c("darkgray",
                                "darkgray",
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
                  ylim = c(1.3, 6)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black")) +
  
  # add text
  annotate(geom = "text",
           label = "Treatment",
           x = 2.3, 
           y = 6.35,
           color = "darkgreen",
           alpha = 0.5,
           size = 5) +
  
  annotate(geom = "text",
           label = "Intrinsic",
           x = 2.35, 
           y = 2.5,
           color = "darkgray",
           alpha = 0.75,
           size = 5)

# model 2
# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(data = covariate.draws.long.2,
                      aes(x = value,
                          y = name,
                          fill = name),
                      color = "#333333",
                      alpha = 0.5,
                      scale = 0.7,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  # 95%
  geom_errorbarh(data = all.ci.2,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = var,
                     color = var),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  # 50%
  geom_errorbarh(data = all.ci.2,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = var,
                     color = var),
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  
  # points
  geom_point(data = all.ci.2,
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
                               "darkgreen",
                               "darkgreen",
                               "darkgreen",
                               "darkgreen")) +
  
  scale_color_manual(values = c("darkgray",
                                "darkgray",
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
                  ylim = c(1.3, 6)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black")) +
  
  # add text
  annotate(geom = "text",
           label = "Treatment",
           x = 2.3, 
           y = 6.35,
           color = "darkgreen",
           alpha = 0.5,
           size = 5) +
  
  annotate(geom = "text",
           label = "Intrinsic",
           x = 2.35, 
           y = 2.5,
           color = "darkgray",
           alpha = 0.75,
           size = 5)

# model 3
# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(data = covariate.draws.long.3,
                      aes(x = value,
                          y = name,
                          fill = name),
                      color = "#333333",
                      alpha = 0.5,
                      scale = 0.7,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  # 95%
  geom_errorbarh(data = all.ci.3,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = var,
                     color = var),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  # 50%
  geom_errorbarh(data = all.ci.3,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = var,
                     color = var),
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  
  # points
  geom_point(data = all.ci.3,
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
                               "darkgreen",
                               "darkgreen",
                               "darkgreen",
                               "darkgreen")) +
  
  scale_color_manual(values = c("darkgray",
                                "darkgray",
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
                  ylim = c(1.3, 6)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black")) +
  
  # add text
  annotate(geom = "text",
           label = "Treatment",
           x = 2.3, 
           y = 6.35,
           color = "darkgreen",
           alpha = 0.5,
           size = 5) +
  
  annotate(geom = "text",
           label = "Intrinsic",
           x = 2.35, 
           y = 2.5,
           color = "darkgray",
           alpha = 0.75,
           size = 5)

#_______________________________________________________________________________________________
# 5. All scenarios ----

# bind together
covariate.draws.long.all <- rbind(covariate.draws.long.1,
                                  covariate.draws.long.2,
                                  covariate.draws.long.3)

# variable name order
covariate.draws.long.all$name <- factor(covariate.draws.long.all$name,
                                        levels = c("Retention (year 1)",
                                                   "Retention (year 2)",
                                                   "Piling (year 1)",
                                                   "Piling (year 2)",
                                                   "Body condition index",
                                                   "Body condition index * t"))

all.ci.all <- rbind(all.ci.1,
                    all.ci.2,
                    all.ci.3)

# change variable name
colnames(all.ci.all)[6] <- "name"

# we'll color and shape each one differently
scen.colors <- c("gray", "red", "gold")
scen.shapes <- c(21, 24, 22)

#_______________________________________________________________________________________________

# plot
ggplot(data = covariate.draws.long.all) +
  
  # white background
  theme_bw() +
  
  # facet
  facet_wrap(~ name,
             ncol = 2,
             strip.position = "left") +
  
  # KDE
  geom_density_ridges(data = covariate.draws.long.all,
                      aes(x = value,
                          y = rev(scen),
                          fill = scen),
                      color = "#333333",
                      alpha = 0.5,
                      scale = 0.7,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  # 95%
  geom_errorbarh(data = all.ci.all,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = rev(scen),
                     color = scen),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  # 50%
  geom_errorbarh(data = all.ci.all,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = rev(scen),
                     color = scen),
                 height = 0,
                 linewidth = 2,
                 position = position_nudge(y = -0.1)) +
  
  # points
  geom_point(data = all.ci.all,
             aes(x = med,
                 y = rev(scen),
                 shape = scen,
                 size = scen),
             color = "black",
             fill = "white",
             stroke = 0.8,
             position = position_nudge(y = -0.1)) +
  
  # colors and shapes
  scale_fill_manual(values = scen.colors,
                    labels = c("Scenario 1",
                               "Scenario 2",
                               "Scenario 3")) +
  
  scale_color_manual(values = scen.colors,
                     labels = c("Scenario 1",
                                "Scenario 2",
                                "Scenario 3")) +
    
  scale_shape_manual(values = scen.shapes,
                     labels = c("Scenario 1",
                                "Scenario 2",
                                "Scenario 3")) +
  
  scale_size_manual(values = c(1.5, 1.25, 1.25),
                    labels = c("Scenario 1",
                               "Scenario 2",
                               "Scenario 3")) +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "gray45") +
  
  # axis titles
  xlab("Hazard ratio") +
  ylab("") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3)) +
  
  # coordinates
  coord_cartesian(xlim = c(0.35, 3.25),
                  ylim = c(1.25, 3.25)) +
  
  # theme arguments
  theme(
    
    # facet strip
    strip.background = element_rect(fill = NA,
                                    color = NA),
    strip.text = element_text(hjust = 1,
                              vjust = 1),
    
    # legend
    legend.position = c(0.85, 0.16),
    legend.title = element_blank(),
    
    # panel border and gridlines
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    
    # remove y-axis ticks and labels
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # black axis text
    axis.text = element_text(color = "black")
    
    )

# 598 x 490

#_______________________________________________________________________________________________
# 6. Write HDIs to clipboard ----
#_______________________________________________________________________________________________

write.table(all.ci.all, "clipboard", sep = "\t")
