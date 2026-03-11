# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 08c - Posterior predictive checks - results
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 27 Feb 2026 
# Date completed: 02 Mar 2026 
# Date last modified: 11 Mar 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

discrep.i.all <- readRDS("PPCs/discrep_i_all.rds")
discrep.sum.all <- readRDS("PPCs/discrep_sum_all.rds")
S.df.all <- readRDS("PPCs/S_df_all.rds")

# fates dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

#_______________________________________________________________________________________________
# 3. Bayesian p-value ----
#_______________________________________________________________________________________________
# 3a. Calculation ----
#_______________________________________________________________________________________________

(bayes.p <- discrep.i.all %>%
   
   summarize(scen1 = sum(Dsim.Dobs.1) / n(),
             scen2 = sum(Dsim.Dobs.2) / n(),
             scen3 = sum(Dsim.Dobs.3) / n())
 
)

#_______________________________________________________________________________________________
# 3b. Plot ----

# mirroring Jones et al. (2020)
# https://wildlife.onlinelibrary.wiley.com/doi/10.1002/jwmg.21886

#_______________________________________________________________________________________________

# pivot
discrep.sum.pivot <- discrep.sum.all %>%
  
  # initial longer pivot
  pivot_longer(cols = c(discrep.sim.1.sum:discrep.obs.3.sum)) %>%
  
  # new columns
  mutate(
    
    scenario = factor(substr(name, 13, 13),
                      labels = c("Scenario 1",
                                 "Scenario 2",
                                 "Scenario 3")),
    
    which.discrep = substr(name, 9, 11)
    
  ) %>%
  
  # remove name
  dplyr::select(-name) %>%
  
  # pivot wider
  pivot_wider(names_from = which.discrep)

# plot
ggplot(data = discrep.sum.pivot) +
  
  theme_bw() +
  
  facet_wrap(~ scenario, nrow = 1) +
  
  # 1:1 line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  
  geom_point(aes(x = obs,
                 y = sim),
             alpha = 0.25,
             size = 0.75,
             shape = 21,
             color = "aquamarine4") +
  
  # axes
  coord_cartesian(xlim = c(90, 275),
                  ylim = c(90, 275)) +
  
  xlab("Observed discrepancy") +
  ylab("Simulated discrepancy") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color ="black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0))

# 600 x 250

#_______________________________________________________________________________________________
# 4. Survivorship curves ----
#_______________________________________________________________________________________________
# 4a. Functions ----
#_______________________________________________________________________________________________

# function to incorporate t since monitoring started
add_t <- function (x) {
  
  x$t = 1:nrow(x)
  
  return(x)
  
}

# Kaplan Meier for empirical
kap_mei <- function (
    
  x,
  response = "y.mort.scen1"
  
) {
  
  # loop through all follow-up times
  all.S <- vector(length = max(x$t))
  s.i <- vector(length = max(x$t))
  all.ci.low <-vector(length = max(x$t))
  all.ci.upp <-vector(length = max(x$t))
  
  all.S[1] = 1.0 # initialize at 100%
  s.i[1] = 1.0 # initialize at 100%
  all.ci.low[1] <- 1.0
  all.ci.upp[1] <- 1.0
  
  for (i in 2:max(x$t)) {
    
    # subset
    indivs.i <- x[x$t == i, ]
    
    # calculate survival at interval
    s.i[i] = 1 - (sum(indivs.i[[response]]) / nrow(indivs.i))
    
    # cumulative survival at interval
    all.S[i] = prod(s.i[1:i])
    
    # calculate 90 % confidence interval
    focal.var = ((all.S[i]^2) * (1 - all.S[i])) / nrow(indivs.i)
    
    all.ci.low[i] = all.S[i] - (1.64 * sqrt(focal.var))
    all.ci.upp[i] = all.S[i] + (1.64 * sqrt(focal.var))
    
  }
  
  # df for returning
  km.df <- data.frame(t = 1:max(x$t),
                      S = all.S,
                      ci.low = all.ci.low,
                      ci.upp = all.ci.upp)
  
  return(km.df)
  
}

# prediction envelopes
km_pred <- function (x) {
  
  x.1 <- x %>%
    
    group_by(scen, t) %>%
    
    mutate(med = median(S),
           pe.low = quantile(S, prob = 0.025),
           pe.upp = quantile(S, prob = 0.975)) %>%
    
    slice(1) %>%
    
    ungroup() %>%
    
    dplyr::select(t, med, pe.low, pe.upp, scen)
  
  return(x.1)
  
}

#_______________________________________________________________________________________________
# 4b. Attribute study-week to fates df ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  # join in year
  mutate(year = fates.year$year) %>%
  
  left_join(
    
    day.lookup %>%
      
      dplyr::select(
        
        year,
        study.year.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  )


%>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1),
         study.week.s = (study.week - mean(study.week)) / sd(study.week),
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________________________
# 4c. Apply functions ----
#_______________________________________________________________________________________________

# empirical data
fates.1.deploy.split <- split(fates.1, fates.1$deployment.1)
fates.2 <- do.call(rbind, lapply(fates.1.deploy.split, add_t))

# KM curves
km.df.all <- rbind(kap_mei(fates.2, response = "y.mort.scen1") %>% mutate(scen = "1"),
                   kap_mei(fates.2, response = "y.mort.scen2") %>% mutate(scen = "2"),
                   kap_mei(fates.2, response = "y.mort.scen3") %>% mutate(scen = "3"))

# simulated curves
km.ppc <- km_pred(S.df.all)  

#_______________________________________________________________________________________________
# 4d. Plot with empirical curve ----

# using all the curves breaks the vector graphics
# the upper and lower prediction envelopes seem reasonable

# define theme
km_theme <- function () {
  
  theme_bw() +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          strip.text = element_text(hjust = 0),
          strip.background = element_rect(fill = "white"))
  
}

# change "scen" labels
km.ppc$scen.name <- paste0("Scenario ", km.ppc$scen)
km.df.all$scen.name <- paste0("Scenario ", km.df.all$scen)

#_______________________________________________________________________________________________

ggplot() +
  
  km_theme() +
  
  # facet
  facet_wrap(~ scen.name) +
  
  # simulated median
  geom_line(data = km.ppc,
            aes(x = t,
                y = med),
            color = "aquamarine4",
            linewidth = 0.65) +
  
  # simulated prediction envelopes
  geom_ribbon(data = km.ppc,
              aes(x = t,
                  y = med,
                  ymin = pe.low,
                  ymax = pe.upp),
              alpha = 0.35,
              fill = "aquamarine4") +
  
  # empirical curve confidence envelope
  geom_ribbon(data = km.df.all,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.15,
              color = "gray50",
              linewidth = 0.35,
              fill = NA,
              linetype = "dashed") +
  
  # empirical curve
  geom_line(data = km.df.all,
            aes(x = t,
                y = S),
            linewidth = 0.35,
            linetype = "dashed") +
  
  xlab("Weeks") +
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(4.25, 70))

# 600 x 250
