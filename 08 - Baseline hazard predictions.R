# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08 - Baseline hazard predictions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 16 Dec 2025 
# Date last modified: 23 Jan 2026 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in and clear data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- readRDS("Model outputs/model_1.rds")
model.fit.2 <- readRDS("Model outputs/model_2.rds")
model.fit.3 <- readRDS("Model outputs/model_3.rds")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# convert to data.frame
model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

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
# 4. Calculate each hazard prediction ----
#_______________________________________________________________________________________________

# function (accepts model fit)
haz_spline <- function (
    
  x,
  low = 0.025,
  upp = 0.975
  
  ) {
  
  # keep only columns we need
  x.1 <- x %>%
    
    dplyr::select(c("a0[1]":"a0[4]",
                    "w[1, 1]":"w[4, 9]"))
  
  # calculate predictions for each sex-forest type group
  # internal function to calculate predictions
  spline_pred <- function (y, group) {
    
    # define a blank matrix [n.draws, n.prediction points]
    pred.matrix <- matrix(data = NA,
                          nrow = nrow(x.1),
                          ncol = nrow(basis.pred))
    
    # convert coefficients to column vector for multiplication
    w <- matrix(data = c(y[[paste0("w[", group, ", 1]")]],
                         y[[paste0("w[", group, ", 2]")]],
                         y[[paste0("w[", group, ", 3]")]],
                         y[[paste0("w[", group, ", 4]")]],
                         y[[paste0("w[", group, ", 5]")]],
                         y[[paste0("w[", group, ", 6]")]],
                         y[[paste0("w[", group, ", 7]")]],
                         y[[paste0("w[", group, ", 8]")]],
                         y[[paste0("w[", group, ", 9]")]]),
           ncol = 9)

    
    # multiply with 'sweep'
    w.by.b <- sweep(basis.pred, 2, w, `*`)
    
    # sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(w.by.b, 1, sum) 
    
    # add to the intercept (a0) and exponentiate to transform to hazard rate scale
    spline.pred <- exp(as.numeric(y[paste0("a0[", group, "]")]) + w.by.b.sum) 
    
    # and add into the array
    return(spline.pred)
    
  }
  
  # F - SFL
  sft.1.spline <- apply(x.1, 1, spline_pred, group = 1)
  
  # summarize
  sft.1.spline.summary <- data.frame(
    
    med = apply(sft.1.spline, 1, median),
    low = apply(sft.1.spline, 1, quantile, prob = low),
    upp = apply(sft.1.spline, 1, quantile, prob = upp),
    sft = 1
    
  )
  
  # M - SFL
  sft.2.spline <- apply(x.1, 1, spline_pred, group = 2)
  
  # summarize
  sft.2.spline.summary <- data.frame(
    
    med = apply(sft.2.spline, 1, median),
    low = apply(sft.2.spline, 1, quantile, prob = low),
    upp = apply(sft.2.spline, 1, quantile, prob = upp),
    sft = 2
    
  )
  
  # F - XMC
  sft.3.spline <- apply(x.1, 1, spline_pred, group = 3)
  
  # summarize
  sft.3.spline.summary <- data.frame(
    
    med = apply(sft.3.spline, 1, median),
    low = apply(sft.3.spline, 1, quantile, prob = low),
    upp = apply(sft.3.spline, 1, quantile, prob = upp),
    sft = 3
    
  )
  
  # M - XMC
  sft.4.spline <- apply(x.1, 1, spline_pred, group = 4)
  
  # summarize
  sft.4.spline.summary <- data.frame(
    
    med = apply(sft.4.spline, 1, median),
    low = apply(sft.4.spline, 1, quantile, prob = low),
    upp = apply(sft.4.spline, 1, quantile, prob = upp),
    sft = 4
    
  )
  
  # bind together for further manipulation
  all.spline.summary <- rbind(sft.1.spline.summary,
                              sft.2.spline.summary,
                              sft.3.spline.summary,
                              sft.4.spline.summary)
  
  # calculate hazard ratios FIRST on the full posterior predictions
  # SFL is 3 / 1
  # XMC is 4 / 2
  SFL.hr.sex <- sft.3.spline / sft.1.spline
  XMC.hr.sex <- sft.4.spline / sft.2.spline
  
  # make a list
  spline.list <- list(all.spline.summary,
                      SFL.hr.sex,
                      XMC.hr.sex)
  
  # return
  return(spline.list)
  
}

# apply function
spline.summary.mod1 <- haz_spline(model.fit.1)
spline.summary.mod2 <- haz_spline(model.fit.2)
spline.summary.mod3 <- haz_spline(model.fit.3)

#_______________________________________________________________________________________________
# 4. Process predictions for each plot ----
#_______________________________________________________________________________________________
# 4a. Mean hazard across the annual cycle, by forest type ----

# let's take the mean prediction for each sex, and show the sex differences as 
# hazard ratios in the next plot

# SFL is c(1, 3)

#_______________________________________________________________________________________________

# function
spline_ft <- function(x) {
  
  
  spline.summary.ft <- x %>%
    
    mutate(
      
      # t
      t = rep(weeks.pred, times = 4),
      
      # forest type identifier
      ft = ifelse(sft %in% c(1, 3),
                  "SFL",
                  "XMC")) %>%
    
    group_by(t, ft) %>%
    
    summarize(med = mean(med),
              low = mean(low),
              upp = mean(upp)) %>%
    
    ungroup()
  
  return(spline.summary.ft)

}

# apply function
spline.ft.mod1 <- spline_ft(spline.summary.mod1[[1]])
spline.ft.mod2 <- spline_ft(spline.summary.mod2[[1]])
spline.ft.mod3 <- spline_ft(spline.summary.mod3[[1]])

#_______________________________________________________________________________________________
# 4b. Hazard ratios by sex ----

# here we'll calculate the M:F hazard ratio across the annual cycle, by forest type

# F is c(1, 3)

#_______________________________________________________________________________________________

# function
spline_sex <- function(
    
  x,
  low = 0.025,
  upp = 0.975
  
  ) {
  
  # SFL
  spline.summary.sex <- bind_rows(
    
    data.frame(
    
    t = weeks.pred,
    hr.med = apply(x[[2]], 1, median),
    hr.low = apply(x[[2]], 1, quantile, prob = low),
    hr.upp = apply(x[[2]], 1, quantile, prob = upp),
    ft = "SFL"
    
  ),
  
  data.frame(
    
    t = weeks.pred,
    hr.med = apply(x[[3]], 1, median),
    hr.low = apply(x[[3]], 1, quantile, prob = low),
    hr.upp = apply(x[[3]], 1, quantile, prob = upp),
    ft = "XMC"
    
  )
  
)
    
  return(spline.summary.sex)
  
}

# apply function
spline.sex.mod1 <- spline_sex(spline.summary.mod1)
spline.sex.mod2 <- spline_sex(spline.summary.mod2)
spline.sex.mod3 <- spline_sex(spline.summary.mod3)

#_______________________________________________________________________________________________
# 5. Plot predictions ----
#_______________________________________________________________________________________________

# define month cutoffs
month.cutoffs <- data.frame(breaks = c(1, 5.57, 9.86, 
                                       14.29, 18.71, 22.71, 
                                       27.14, 31.43, 35.86, 
                                       40.14, 44.57, 48.86),
                            labels = c("O", "N", "D", 
                                       "J", "F", "M", 
                                       "A", "M", "J", 
                                       "J", "A", "S"))

#_______________________________________________________________________________________________
# 5a. Mean hazard across the annual cycle, by forest type ----
#_______________________________________________________________________________________________

# scenario 1
ggplot(data = 
         
         spline.ft.mod1 %>% mutate(ft = factor(ft,
                                               levels = c("SFL", "XMC"),
                                               labels = c("a", "b")))
       
       ) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = low,
                  ymax = upp,
                  fill = ft),
              alpha = 0.25) +
  
  # median
  geom_line(aes(x = t,
                y = med,
                color = ft),
            linewidth = 1.5) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)) +
  
  # axis labels
  ylab("Baseline hazard") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # colors
  scale_color_manual(values = c("#003300", "#669900")) +
  scale_fill_manual(values = c("#003300", "#669900"))

# 342 x 589

# scenario 2
ggplot(data = 
         
         spline.ft.mod2 %>% mutate(ft = factor(ft,
                                               levels = c("SFL", "XMC"),
                                               labels = c("a", "b")))
       
) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = low,
                  ymax = upp,
                  fill = ft),
              alpha = 0.25) +
  
  # median
  geom_line(aes(x = t,
                y = med,
                color = ft),
            linewidth = 1.5) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)) +
  
  # axis labels
  ylab("Baseline hazard") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # colors
  scale_color_manual(values = c("#003300", "#669900")) +
  scale_fill_manual(values = c("#003300", "#669900"))

# 342 x 589

# scenario 3
ggplot(data = 
         
         spline.ft.mod3 %>% mutate(ft = factor(ft,
                                               levels = c("SFL", "XMC"),
                                               labels = c("a", "b")))
       
) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = low,
                  ymax = upp,
                  fill = ft),
              alpha = 0.25) +
  
  # median
  geom_line(aes(x = t,
                y = med,
                color = ft),
            linewidth = 1.5) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)) +
  
  # axis labels
  ylab("Baseline hazard") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # colors
  scale_color_manual(values = c("#003300", "#669900")) +
  scale_fill_manual(values = c("#003300", "#669900"))

# 342 x 589

#_______________________________________________________________________________________________
# 5b. M:F hazard ratio across the annual cycle, by forest type ----
#_______________________________________________________________________________________________

# scenario 1
ggplot(data = 
         
         spline.sex.mod1 %>% mutate(
           
           # forest type factor labels
           ft = factor(ft,
                       levels = c("SFL", "XMC"),
                       labels = c("a", "b")),
           
           # conditional colors
           Color = factor(ifelse(hr.med > 1,
                                 "M",
                                 "F"),
                          levels = c("M", "F"),
                          labels = c("male HR > 1",
                                     "male HR < 1"))
           
         )
       
) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # horizontal line at HR == 1
  geom_hline(yintercept = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = hr.low,
                  ymax = hr.upp),
              alpha = 0.2) +
  
  # medians
  # underlying male line for continuity
  geom_line(aes(x = t,
                y = hr.med),
            linewidth = 1.5,
            color = "gray25") +
  
  # here' well split datasets based on whether the median is above or below 1
  # males
  geom_line(
    
    data = spline.sex.mod1 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med >= 1),
    
    aes(x = t,
        y = hr.med),
        linewidth = 1.5,
    color = "gray25") +
  
  # females (first segment)
  geom_line(
    
    data = spline.sex.mod1 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t < 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # females (second segment)
  geom_line(
    
    data = spline.sex.mod1 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t > 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = c(0.7, 0.4),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # point size in legend
  guides(color = guide_legend(override.aes = list(size = 3))) +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  # axis labels
  ylab("M:F hazard ratio") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7))

# 342 x 589

# scenario 2
ggplot(data = 
         
         spline.sex.mod2 %>% mutate(
           
           # forest type factor labels
           ft = factor(ft,
                       levels = c("SFL", "XMC"),
                       labels = c("a", "b")),
           
           # conditional colors
           Color = factor(ifelse(hr.med > 1,
                                 "M",
                                 "F"),
                          levels = c("M", "F"),
                          labels = c("male HR > 1",
                                     "male HR < 1"))
           
         )
       
) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # horizontal line at HR == 1
  geom_hline(yintercept = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = hr.low,
                  ymax = hr.upp),
              alpha = 0.2) +
  
  # medians
  # underlying male line for continuity
  geom_line(aes(x = t,
                y = hr.med),
            linewidth = 1.5,
            color = "gray25") +
  
  # here' well split datasets based on whether the median is above or below 1
  # males (first segment)
  geom_line(
    
    data = spline.sex.mod2 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med >= 1 & t < 10),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "gray25") +
  
  # males (second segment)
  geom_line(
    
    data = spline.sex.mod2 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med >= 1 & t > 10 & t < 40),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "gray25") +
  
  # females (first segment)
  geom_line(
    
    data = spline.sex.mod2 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t < 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # females (second segment)
  geom_line(
    
    data = spline.sex.mod2 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t > 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = c(0.7, 0.4),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # point size in legend
  guides(color = guide_legend(override.aes = list(size = 3))) +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  # axis labels
  ylab("M:F hazard ratio") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7))

# 342 x 589

# scenario 3
ggplot(data = 
         
         spline.sex.mod3 %>% mutate(
           
           # forest type factor labels
           ft = factor(ft,
                       levels = c("SFL", "XMC"),
                       labels = c("a", "b")),
           
           # conditional colors
           Color = factor(ifelse(hr.med > 1,
                                 "M",
                                 "F"),
                          levels = c("M", "F"),
                          labels = c("male HR > 1",
                                     "male HR < 1"))
           
         )
       
) +
  
  theme_bw() +
  
  facet_wrap(~ ft,
             ncol = 1) +
  
  # horizontal line at HR == 1
  geom_hline(yintercept = 1) +
  
  # LIGHT vertical lines for minor month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "gray",
             linetype = "dashed",
             alpha = 0.5) +
  
  # vertical lines for major month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(2, 4, 6, 8, 10, 12)],
             color = "gray") +
  
  # credible interval
  geom_ribbon(aes(x = t,
                  ymin = hr.low,
                  ymax = hr.upp),
              alpha = 0.2) +
  
  # medians
  # underlying male line for continuity
  geom_line(aes(x = t,
                y = hr.med),
            linewidth = 1.5,
            color = "gray25") +
  
  # here' well split datasets based on whether the median is above or below 1
  # males (first segment)
  geom_line(
    
    data = spline.sex.mod3 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med >= 1 & t < 15),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "gray25") +
  
  # males (second segment)
  geom_line(
    
    data = spline.sex.mod3 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med >= 1 & t > 15),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "gray25") +
  
  # females (first segment)
  geom_line(
    
    data = spline.sex.mod3 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t < 10),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # females (second segment)
  geom_line(
    
    data = spline.sex.mod3 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t > 15 & t < 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # females (third segment)
  geom_line(
    
    data = spline.sex.mod3 %>% mutate(
      
      # forest type factor labels
      ft = factor(ft,
                  levels = c("SFL", "XMC"),
                  labels = c("a", "b"))
      
    ) %>%
      
      # filter
      filter(hr.med <= 1 & t > 30),
    
    aes(x = t,
        y = hr.med),
    linewidth = 1.5,
    color = "#FF3300") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0.25,
                                   face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = c(0.7, 0.4),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black"),
        
        # text inside panels
        strip.text = element_text(hjust = -0,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  # point size in legend
  guides(color = guide_legend(override.aes = list(size = 3))) +
  
  # axis labels
  scale_x_continuous(breaks = month.cutoffs$breaks + 2,
                     labels = month.cutoffs$labels) +
  
  # axis labels
  ylab("M:F hazard ratio") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7))

# 342 x 589

