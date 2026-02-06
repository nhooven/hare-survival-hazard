# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08 - Effect size plots
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2025 
# Date completed: 01 Dec 2025
# Date last modified: 04 Feb 2026
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
model.fit.1 <- readRDS("Model outputs/model_1.rds")
model.fit.2 <- readRDS("Model outputs/model_2.rds")
model.fit.3 <- readRDS("Model outputs/model_3.rds")

# convert to data.frame
model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

#_______________________________________________________________________________________________
# 3. Visualize covariate effects and evaluate strength of evidence ----
#_______________________________________________________________________________________________
# 3a. Select hazard ratios and pivot_longer ----

# function
hr_pivot <- function (x, scen = "1") {
  
  x.1 <- x %>%
    
    dplyr::select(hr_bci,
                  hr_bci_study_week,
                  hr_dm,
                  hr_open,
                  hr_ret_total1,
                  hr_ret_total2,
                  hr_pil_total1,
                  hr_pil_total2) %>%
    
    pivot_longer(cols = hr_bci:hr_pil_total2) %>%
    
    # reorder factor
    mutate(name = factor(name,
                         levels = c("hr_bci",
                                    "hr_bci_study_week",
                                    "hr_dm",
                                    "hr_open",
                                    "hr_ret_total1",
                                    "hr_ret_total2",
                                    "hr_pil_total1",
                                    "hr_pil_total2"),
                         labels = c("BCI",
                                    "BCI * t", 
                                    "DM",
                                    "OPEN",
                                    "RET (yr 1)",
                                    "RET (yr 2)",
                                    "PIL (yr 1)",
                                    "PIL (yr 2)")),
           scen = scen)
  
  return(x.1)
  
}

#_______________________________________________________________________________________________

# use function
long.draws.1 <- hr_pivot(model.fit.1, "1")
long.draws.2 <- hr_pivot(model.fit.2, "2")
long.draws.3 <- hr_pivot(model.fit.3, "3")

#_______________________________________________________________________________________________
# 3b. Extract median and credible intervals ----

# function
extract_med_ci <- function (x, ci.1 = 0.50, ci.2 = 0.95) {
  
  
  x.1 <- x %>% 
    
    group_by(name) %>%
    
    summarize(
      
      med = median(value),
      sd = sd(value),
      lo.1 = as.numeric(hdi(value, ci = ci.1)[2]),
      up.1 = as.numeric(hdi(value, ci = ci.1)[3]),
      lo.2 = as.numeric(hdi(value, ci = ci.2)[2]),
      up.2 = as.numeric(hdi(value, ci = ci.2)[3])
      
    ) %>%
    
    ungroup()
  
  return (x.1)
  
}

#_______________________________________________________________________________________________

# use function
all.ci.1 <- extract_med_ci(long.draws.1)
all.ci.2 <- extract_med_ci(long.draws.2)
all.ci.3 <- extract_med_ci(long.draws.3)

# add scen
all.ci.1$scen = "1"
all.ci.2$scen = "2"
all.ci.3$scen = "3"

#_______________________________________________________________________________________________
# 4. Half-eye plots - by scenario ----
#_______________________________________________________________________________________________
# 4a. Function ----
#_______________________________________________________________________________________________

plot_by_scen <- function (long.draws, ci) {
  
  # add index for type of variable
  long.draws.forPlot <- long.draws %>%
    
    mutate(
      
     vartype = case_when(
       
       name %in% c("BCI", "BCI * t") ~ "intrinsic",
       name %in% c("DM", "OPEN") ~ "landscape",
       name %in% c("RET (yr 1)",
                   "RET (yr 2)",
                   "PIL (yr 1)",
                   "PIL (yr 2)") ~ "treatment"
       
     )
      
  )
  
  ci.forPlot <- ci %>%
    
    mutate(
      
      vartype = case_when(
        
        name %in% c("BCI", "BCI * t") ~ "intrinsic",
        name %in% c("DM", "OPEN") ~ "landscape",
        name %in% c("RET (yr 1)",
                    "RET (yr 2)",
                    "PIL (yr 1)",
                    "PIL (yr 2)") ~ "treatment"
        
      )
      
    )
  
out.plot <- ggplot() +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  # KDE
  geom_density_ridges(data = long.draws.forPlot,
                      aes(x = value,
                          y = name,
                          fill = vartype),
                      color = NA,
                      alpha = 0.5,
                      scale = 0.5,
                      rel_min_height = 0.006) +     # remove tails
  
  # credible intervals
  geom_errorbarh(data = ci.forPlot,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = vartype),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.75,
                 position = position_nudge(y = -0.1)) +
  
  geom_errorbarh(data = ci.forPlot,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = vartype),
                 height = 0,
                 linewidth = 1.75,
                 position = position_nudge(y = -0.1)) +
  
  # fill and color
  scale_fill_manual(values = c("gray60", "green4", "#FF3300")) +
  scale_color_manual(values = c("gray40", "green4", "#FF3300")) +
  
  # axis title
  xlab("Hazard ratio") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3)) +
  
  # y-axis scale
  scale_y_discrete(limits = rev) +
  
  # coordinates
  coord_cartesian(xlim = c(0.45, 2.95),
                  ylim = c(1.3, 8)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        strip.background = element_blank())
    
  return(out.plot)
  
}

#_______________________________________________________________________________________________
# 4b. Use ----

# 379 x 420

#_______________________________________________________________________________________________

(hr.scen1 <- plot_by_scen(long.draws.1, all.ci.1))
(hr.scen2 <- plot_by_scen(long.draws.2, all.ci.2))
(hr.scen3 <- plot_by_scen(long.draws.3, all.ci.3))

#_______________________________________________________________________________________________
# 5. Credible interval plots - All scenarios ----

# this will only be the CIs
all.ci.all <- rbind(all.ci.1,
                    all.ci.2,
                    all.ci.3)

# we'll color and fill each one differently
scen.colors <- c("black", "red", "gold")

#_______________________________________________________________________________________________

# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  # credible intervals
  geom_errorbarh(data = all.ci.all,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = scen),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  geom_errorbarh(data = all.ci.all,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = scen),
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  # fill and color
  scale_fill_manual("Scenario", values = scen.colors) +
  scale_color_manual("Scenario", values = scen.colors) +
  
  # axis title
  xlab("Hazard ratio") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3)) +
  
  # y-axis scale
  scale_y_discrete(limits = rev) +
  
  # coordinates
  coord_cartesian(xlim = c(0.45, 2.55)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.title = element_text(),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        strip.background = element_blank())

# 379 x 420

#_______________________________________________________________________________________________
# 6. Write HDIs to clipboard ----
#_______________________________________________________________________________________________

write.table(all.ci.all, "clipboard", sep = "\t")
