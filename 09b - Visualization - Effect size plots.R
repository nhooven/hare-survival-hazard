# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 09b - Effect size plots
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 19 Nov 2025 
# Date completed: 01 Dec 2025
# Date last modified: 01 May 2026
# R version: 4.4.3

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(ggridges)        # ridgeline plot
library(bayestestR)      # highest density posterior intervals
library(cowplot)

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

# model samples
model.fit.1 <- readRDS("models/model_1.rds")
model.fit.2 <- readRDS("models/model_2.rds")
model.fit.3 <- readRDS("models/model_3.rds")

# convert to data.frame
model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

# select betas only
b.names <- c("b_bci",
             "b_dm",
             "b_open",
             "b_post1",
             "b_post2",
             "b_ret_post1",
             "b_ret_post2",
             "b_pil_post1",
             "b_pil_post2")

model.fit.1 <- model.fit.1 %>% dplyr::select(all_of(b.names))
model.fit.2 <- model.fit.2 %>% dplyr::select(all_of(b.names))
model.fit.3 <- model.fit.3 %>% dplyr::select(all_of(b.names))

#_______________________________________________________________________________
# 3. Calculate correct hazard ratios ----

# this is easy for bci/dm/open, but will need to be year-specific for the treatments
# remember that the "total" effects are with respect to unthinned/pre-treatment
# we can subtract out the year effect to make an interesting plot
# we'll also pivot for summary and plotting here

#_______________________________________________________________________________
# 3a. Function ----
#_______________________________________________________________________________

calc_hr <- function (x, scen = "1") {
  
  x.1 <- x %>%
    
    mutate(
      
      # easy ones
      hr.bci = exp(b_bci),
      hr.dm = exp(b_dm),
      hr.open = exp(b_open),
      
      # years alone
      hr.post1 = exp(b_post1),
      hr.post2 = exp(b_post2),
      
      # total treatment effects
      hr.ret.post1.total = exp(b_post1 + b_ret_post1),
      hr.ret.post2.total = exp(b_post2 + b_ret_post2),
      hr.pil.post1.total = exp(b_post1 + b_pil_post1),
      hr.pil.post2.total = exp(b_post2 + b_pil_post2),
      
      # year comparison effects
      hr.ret.post1.comp = exp(b_ret_post1 - b_post1),
      hr.ret.post2.comp = exp(b_ret_post2 - b_post2),
      hr.pil.post1.comp = exp(b_pil_post1 - b_post1),
      hr.pil.post2.comp = exp(b_pil_post2 - b_post2)
      
    ) %>%
    
    # pivot
    pivot_longer(hr.bci:hr.pil.post2.comp) %>%
    
    # reorder factor
    mutate(name = factor(name,
                         levels = c("hr.bci",
                                    "hr.dm",
                                    "hr.open",
                                    "hr.post1",
                                    "hr.post2",
                                    "hr.ret.post1.total",
                                    "hr.ret.post2.total",
                                    "hr.pil.post1.total",
                                    "hr.pil.post2.total",
                                    "hr.ret.post1.comp",
                                    "hr.ret.post2.comp",
                                    "hr.pil.post1.comp",
                                    "hr.pil.post2.comp")),
           scen = scen
           
           ) %>%
    
    # keep only relevant columns
    dplyr::select(name, value, scen)
          
  return(x.1)
  
}

# use function
long.draws.1 <- calc_hr(model.fit.1, "1")
long.draws.2 <- calc_hr(model.fit.2, "2")
long.draws.3 <- calc_hr(model.fit.3, "3")

#_______________________________________________________________________________
# 3b. Extract median and credible intervals ----

# function
extract_med_ci <- function (x, ci.1 = 0.50, ci.2 = 0.90) {
  
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

#_______________________________________________________________________________

# use function
all.ci.1 <- extract_med_ci(long.draws.1)
all.ci.2 <- extract_med_ci(long.draws.2)
all.ci.3 <- extract_med_ci(long.draws.3)

# add scen
all.ci.1$scen = "1"
all.ci.2$scen = "2"
all.ci.3$scen = "3"

#_______________________________________________________________________________
# 4. Half-eye plots - by scenario ----

# I envision 3 panels:
# top: bci, dm, open HRs
# middle: year-to-year total effects, split by treatment
# bottom: year-to-year comparative effects

# var lists
panel.1.vars <- all.ci.1$name[1:3]
panel.2.vars <- all.ci.1$name[4:9]
panel.3.vars <- all.ci.1$name[10:13]

#_______________________________________________________________________________
# 4a. Function ----
#_______________________________________________________________________________

plot_by_scen <- function (long.draws, ci) {
  
  # plotlist  
  out.plot <- list()  
  
  # panel 1 - basic HRs
  # prep dfs
  long.draws.panel.1 <- long.draws %>% 
    
    filter(name %in% panel.1.vars) %>%
    
    mutate(name = factor(name,
                         labels = c("BCI", 
                                    "DM",
                                    "OPEN")))
  
  ci.panel.1 <- ci %>% 
    
    filter(name %in% panel.1.vars) %>%
    
    mutate(name = factor(name,
                         labels = c("BCI", 
                                    "DM",
                                    "OPEN")))
  
  # plot
  out.plot[[1]] <- ggplot() +
    
    # white background
    theme_bw() +
    
    # vertical line at 1
    geom_vline(xintercept = 1,
               linetype = "dashed",
               color = "black") +
    
    # KDE
    geom_density_ridges(data = long.draws.panel.1,
                        aes(x = value,
                            y = name,
                            fill = name),
                        color = NA,
                        alpha = 0.5,
                        scale = 0.5,
                        rel_min_height = 0.006) +     # remove tails
    
    # credible intervals
    geom_errorbarh(data = ci.panel.1,
                   aes(xmin = lo.2,
                       xmax = up.2,
                       y = name,
                       color = name),
                   alpha = 0.40,
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    geom_errorbarh(data = ci.panel.1,
                   aes(xmin = lo.1,
                       xmax = up.1,
                       y = name,
                       color = name),
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    # x-axis scale
    scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
    
    # y-axis scale
    scale_y_discrete(limits = rev) +
    
    # coordinates
    coord_cartesian(xlim = c(0.27, 4.3),
                    ylim = c(1.3, 3)) +
    
    # remove gridlines, remove legend
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          plot.margin = margin(0.1, 0.1, 0.1, 0.83, unit = "cm")) +
    
    # colors
    scale_color_manual(values = c("gray45", "olivedrab", "olivedrab")) +
    scale_fill_manual(values = c("gray45", "olivedrab", "olivedrab")) +
    
    # annotation
    annotate("text", x = 4.3, y = 3.4, label = "a")
  
  # panel 2 - total effect HRs
  # prep dfs
  long.draws.panel.2 <- long.draws %>% 
    
    filter(name %in% panel.2.vars) %>%
    
    # name order and labels
    mutate(name = factor(name,
                         levels = c("hr.post1",
                                    "hr.ret.post1.total",
                                    "hr.pil.post1.total",
                                    "hr.post2",
                                    "hr.ret.post2.total",
                                    "hr.pil.post2.total"),
                         labels = c("YR 1",
                                    "YR 1 (RET)",
                                    "YR 1 (PIL)",
                                    "YR 2",
                                    "YR 2 (RET)",
                                    "YR 2 (PIL)")))
  
  ci.panel.2 <- ci %>% 
    
    filter(name %in% panel.2.vars) %>%
    
    # name order and labels
    mutate(name = factor(name,
                         levels = c("hr.post1",
                                    "hr.ret.post1.total",
                                    "hr.pil.post1.total",
                                    "hr.post2",
                                    "hr.ret.post2.total",
                                    "hr.pil.post2.total"),
                         labels = c("YR 1",
                                    "YR 1 (RET)",
                                    "YR 1 (PIL)",
                                    "YR 2",
                                    "YR 2 (RET)",
                                    "YR 2 (PIL)")))
  
  # plot
  out.plot[[2]] <- ggplot() +
    
    # white background
    theme_bw() +
    
    # vertical line at 1
    geom_vline(xintercept = 1,
               linetype = "dashed",
               color = "black") +
    
    # KDE
    geom_density_ridges(data = long.draws.panel.2,
                        aes(x = value,
                            y = name,
                            fill = name),
                        color = NA,
                        alpha = 0.5,
                        scale = 0.5,
                        rel_min_height = 0.006) +     # remove tails
    
    # credible intervals
    geom_errorbarh(data = ci.panel.2,
                   aes(xmin = lo.2,
                       xmax = up.2,
                       y = name,
                       color = name),
                   alpha = 0.40,
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    geom_errorbarh(data = ci.panel.2,
                   aes(xmin = lo.1,
                       xmax = up.1,
                       y = name,
                       color = name),
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    # colors
    scale_color_manual(values = c("gray45", "firebrick", "firebrick",
                                  "gray45", "firebrick", "firebrick")) +
    
    scale_fill_manual(values = c("gray45", "firebrick", "firebrick",
                                 "gray45", "firebrick", "firebrick")) +
    
    # x-axis scale
    scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
    
    # y-axis scale
    scale_y_discrete(limits = rev) +
    
    # coordinates
    coord_cartesian(xlim = c(0.27, 4.3),
                    ylim = c(1.1, 6)) +
    
    # remove gridlines, remove legend
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")) +
    
    # annotation
    annotate("text", x = 4.3, y = 6.3, label = "b")
  
  # panel 3 - comparison effects
  # prep dfs
  long.draws.panel.3 <- long.draws %>% 
    
    filter(name %in% panel.3.vars) %>%
    
    # name order and labels
    mutate(name = factor(name,
                         levels = c("hr.ret.post1.comp",
                                    "hr.pil.post1.comp",
                                    "hr.ret.post2.comp",
                                    "hr.pil.post2.comp"),
                         labels = c("RET / YR 1",
                                    "PIL / YR 1",
                                    "RET / YR 2",
                                    "PIL / YR 2")))
  
  ci.panel.3 <- ci %>% 
    
    filter(name %in% panel.3.vars) %>%
    
    # name order and labels
    mutate(name = factor(name,
                         levels = c("hr.ret.post1.comp",
                                    "hr.pil.post1.comp",
                                    "hr.ret.post2.comp",
                                    "hr.pil.post2.comp"),
                         labels = c("RET / YR 1",
                                    "PIL / YR 1",
                                    "RET / YR 2",
                                    "PIL / YR 2")))
  
  out.plot[[3]] <- ggplot() +
    
    # white background
    theme_bw() +
    
    # vertical line at 1
    geom_vline(xintercept = 1,
               linetype = "dashed",
               color = "black") +
    
    # credible intervals
    geom_errorbarh(data = ci.panel.3,
                   aes(xmin = lo.2,
                       xmax = up.2,
                       y = name),
                   color = "firebrick",
                   alpha = 0.40,
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    geom_errorbarh(data = ci.panel.3,
                   aes(xmin = lo.1,
                       xmax = up.1,
                       y = name),
                   color = "firebrick",
                   height = 0,
                   linewidth = 1.75,
                   position = position_nudge(y = -0.1)) +
    
    # axis title
    xlab("Hazard ratio") +
    
    # x-axis scale
    scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
    
    # y-axis scale
    scale_y_discrete(limits = rev) +
    
    # coordinates
    coord_cartesian(xlim = c(0.27, 4.3),
                    ylim = c(1.1, 3.7)) +
    
    # remove gridlines, remove legend
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          axis.text = element_text(color = "black"),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          plot.margin = margin(0.1, 0.1, 0.1, 0.13, unit = "cm")) +
    
    # annotation
    annotate("text", x = 4.3, y = 4, label = "c")
  
  # return full plot
  full.plot <- plot_grid(
    
    plotlist = out.plot, 
    nrow = 3,
    rel_heights = c(0.5, 1.0, 0.5)
    
    )
    
  return(full.plot)
  
}

#_______________________________________________________________________________
# 4b. Use ----

# 341 x 502

#_______________________________________________________________________________

(hr.scen1 <- plot_by_scen(long.draws.1, all.ci.1))
(hr.scen2 <- plot_by_scen(long.draws.2, all.ci.2))
(hr.scen3 <- plot_by_scen(long.draws.3, all.ci.3))

#_______________________________________________________________________________
# 5. Credible interval plots - All scenarios ----

# this will only be the CIs
all.ci.all <- rbind(all.ci.1,
                    all.ci.2,
                    all.ci.3)

# we'll color and fill each one differently
scen.colors <- c("black", "red", "gold")

#_______________________________________________________________________________

# plotlist  
scen.plot <- list()  

# panel 1 - basic HRs
# prep dfs
ci.panel.1 <- all.ci.all %>% 
  
  filter(name %in% panel.1.vars) %>%
  
  mutate(name = factor(name,
                       labels = c("BCI", 
                                  "DM",
                                  "OPEN")))

# plot
scen.plot[[1]] <- ggplot() +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  # credible intervals
  geom_errorbarh(data = ci.panel.1,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = scen),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  geom_errorbarh(data = ci.panel.1,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = scen),
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  
  # y-axis scale
  scale_y_discrete(limits = rev) +
  
  # coordinates
  coord_cartesian(xlim = c(0.27, 4.3),
                  ylim = c(1.1, 2.8)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.83, unit = "cm")) +
  
  scale_color_manual(values = scen.colors) +
  
  # annotation
  annotate("text", x = 4.3, y = 3.2, label = "a")

# panel 2 - total effect HRs
# prep dfs
ci.panel.2 <- all.ci.all %>% 
  
  filter(name %in% panel.2.vars) %>%
  
  # name order and labels
  mutate(name = factor(name,
                       levels = c("hr.post1",
                                  "hr.ret.post1.total",
                                  "hr.pil.post1.total",
                                  "hr.post2",
                                  "hr.ret.post2.total",
                                  "hr.pil.post2.total"),
                       labels = c("YR 1",
                                  "YR 1 (RET)",
                                  "YR 1 (PIL)",
                                  "YR 2",
                                  "YR 2 (RET)",
                                  "YR 2 (PIL)")))

# plot
scen.plot[[2]] <- ggplot() +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  # credible intervals
  geom_errorbarh(data = ci.panel.2,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = scen),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  geom_errorbarh(data = ci.panel.2,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = scen),
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  
  # y-axis scale
  scale_y_discrete(limits = rev) +
  
  # coordinates
  coord_cartesian(xlim = c(0.27, 4.3),
                  ylim = c(1.1, 5.9)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = c(0.75, 0.55),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")) +
  
  scale_color_manual(values = scen.colors) +
  
  labs(color = "Scenario") +
  
  # annotation
  annotate("text", x = 4.3, y = 6.3, label = "b")

# panel 3 - comparison effects
# prep dfs
ci.panel.3 <- all.ci.all %>% 
  
  filter(name %in% panel.3.vars) %>%
  
  # name order and labels
  mutate(name = factor(name,
                       levels = c("hr.ret.post1.comp",
                                  "hr.pil.post1.comp",
                                  "hr.ret.post2.comp",
                                  "hr.pil.post2.comp"),
                       labels = c("RET / YR 1",
                                  "PIL / YR 1",
                                  "RET / YR 2",
                                  "PIL / YR 2")))

scen.plot[[3]] <- ggplot() +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  # credible intervals
  geom_errorbarh(data = ci.panel.3,
                 aes(xmin = lo.2,
                     xmax = up.2,
                     y = name,
                     color = scen),
                 alpha = 0.40,
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  geom_errorbarh(data = ci.panel.3,
                 aes(xmin = lo.1,
                     xmax = up.1,
                     y = name,
                     color = scen),
                 height = 0,
                 linewidth = 1.75,
                 position = ggstance::position_dodgev(height = 0.5)) +
  
  # axis title
  xlab("Hazard ratio") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  
  # y-axis scale
  scale_y_discrete(limits = rev) +
  
  # coordinates
  coord_cartesian(xlim = c(0.27, 4.3),
                  ylim = c(1.1, 3.7)) +
  
  # remove gridlines, remove legend
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.13, unit = "cm")) +
  
  scale_color_manual(values = scen.colors) +
  
  # annotation
  annotate("text", x = 4.3, y = 4, label = "c")

# return full plot
plot_grid(
  
  plotlist = scen.plot, 
  nrow = 3,
  rel_heights = c(0.25, 0.55, 0.4)
  
)

# 358 x 502

#_______________________________________________________________________________
# 6. Write HDIs to clipboard ----
#_______________________________________________________________________________

write.table(all.ci.all, "clipboard", sep = "\t")
