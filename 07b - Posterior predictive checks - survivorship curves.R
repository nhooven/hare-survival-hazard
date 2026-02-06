# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06b - Posterior predictive checks - survivorship curves
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 06 Feb 2026 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# Simulated lifetimes
# Here we'll simulate lifetimes from a subset of posterior draws, create Kaplan-Meier
# survivorship curves for each "population", then compare to the empirical K-M curve
# I'll just need to snap every individual to the same starting place

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)
library(tictoc)          # timing

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- readRDS("Model outputs/model_1.rds")
model.fit.2 <- readRDS("Model outputs/model_2.rds")
model.fit.3 <- readRDS("Model outputs/model_3.rds")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________________________
# 3a. Convert mcmc.lists to dfs ----
#_______________________________________________________________________________________________

model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

#_______________________________________________________________________________________________
# 3b. Attribute study-week to fates df ----
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
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1),
         study.week.s = (study.week - mean(study.week)) / sd(study.week),
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________________________
# 3c. Spline for prediction ----
#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# weeks to predict on
weeks.pred <- 1:52

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________________________
# 4. Split dataset by deployment ----
#_______________________________________________________________________________________________

# split by deployment
fates.1.deploy.split <- split(fates.1, fates.1$deployment.1)

#_______________________________________________________________________________________________
# 5. Prepare dataset for Kaplan-Meier curves ----
#_______________________________________________________________________________________________
# 5a. Prepare empirical data ----

# Here we'll keep everything the same, but include a "t" variable denoting the week from the beginning of 
# monitoring per deployment

#_______________________________________________________________________________________________

# function to incorporate t since monitoring started
add_t <- function (x) {
  
  x$t = 1:nrow(x)
  
  return(x)
  
}

# apply
fates.2 <- do.call(rbind, lapply(fates.1.deploy.split, add_t))

#_______________________________________________________________________________________________
# 5b. Write function ----

# Ideally we can also use this within the simulation procedure
# we'll need to split data based on t since deployment

#_______________________________________________________________________________________________

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

# apply function
# model 1
km.df.1 <- kap_mei(fates.2,
                   response = "y.mort.scen1")

# model 2
km.df.2 <- kap_mei(fates.2,
                   response = "y.mort.scen2")

# model 3
km.df.3 <- kap_mei(fates.2,
                   response = "y.mort.scen3")

# bind together
km.df.all <- rbind(km.df.1 %>% mutate(scen = "1"),
                   km.df.2 %>% mutate(scen = "2"),
                   km.df.3 %>% mutate(scen = "3"))

# plot each
ggplot() +
  
  theme_bw() +
  
  # scen1
  geom_ribbon(data = km.df.1,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "gray") +
  
  geom_line(data = km.df.1,
            aes(x = t,
                y = S)) +
  
  # scen2
  geom_ribbon(data = km.df.2,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "red") +
  
  geom_line(data = km.df.2,
            aes(x = t,
                y = S),
            color = "darkred") +
  
  # scen3
  geom_ribbon(data = km.df.3,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "gold") +
  
  geom_line(data = km.df.3,
            aes(x = t,
                y = S),
            color = "orange") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1))

#_______________________________________________________________________________________________
# 6. Simulate lifetimes ----

# Here we'll loop through each iteration (take a subsample first),
# loop through each deployment (n = 383),
# applying a calc intensity function on an 150-week dataset,
# take a Poisson draw and truncate the dataset after the first event
# then apply the K-M curve function

# the output will be a collection of K-M curves

#_______________________________________________________________________________________________
# 6a. Helper functions ----
#_______________________________________________________________________________________________

# expand dataset for max.t weeks
expand_time <- function (x, max.t = 150) {
  
  # x is an individual deployment
  # keep only one row
  x.1 <- x %>%
    
    dplyr::select(deployment.1, site, sex_forest, p.dm.s, p.open.s, BCI.s, study.week) %>%
    
    slice(1) %>%
    
    # named Site
    mutate(
      
      Site = case_when(
        
        site == 1 ~ "1A",
        site == 2 ~ "1B",
        site == 3 ~ "1C",
        site == 4 ~ "2A",
        site == 5 ~ "2B",
        site == 6 ~ "2C",
        site == 7 ~ "3A",
        site == 8 ~ "3B",
        site == 9 ~ "3C",
        site == 10 ~ "4A",
        site == 11 ~ "4B",
        site == 12 ~ "4C"
        
      )
      
    )
  
  # correct study.week (max is 160)
  # we should replace anything above 160 with 160
  t = 1:max.t
  
  # actual study.week
  study.week.real = seq(x.1$study.week, x.1$study.week + (160 - x.1$study.week))
  
  # add concat the max to fill in the rest of the max.t
  # only if necessary
  if (max.t - length(study.week.real) > 0) {
    
    study.week <- c(study.week.real, rep(160, times = max.t - length(study.week.real)))
    
  } else {
    
    # truncate if needed
    study.week <- study.week.real[1:150]
    
  }
  
  # and standardize
  study.week.s <- (study.week - mean(fates.2$study.week)) / sd(fates.2$study.week)
  
  # year.week for starting study.week (take only the first one)
  year.week.start <- day.lookup$year.week[day.lookup$study.week == study.week[1]][1]
  
  # all correct year.weeks
  # this is probably stupid code
  week <- c(seq(year.week.start, 52),
            rep(1:52, times = 5))
  
  # keep only the first max.t  
  week <- week[1:max.t]
  
  # bind into df
  suppressMessages(
    
    x.df <- data.frame(
      
      deploy = x.1$deployment.1,
      Site = x.1$Site,
      sex_forest = x.1$sex_forest,
      BCI.s = x.1$BCI.s,
      p.dm.s = x.1$p.dm.s,
      p.open.s = x.1$p.open.s,
      study.week = study.week,
      study.week.s = study.week.s,
      week = week,
      t = t
      
    ) %>%
      
      # and add treatment variables
      mutate(
        
        # post 1
        post1 = case_when(
          
          # 2 and 3
          # pre-treatment
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week < 53 ~ 0,
          
          # post-treatment 1
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week >= 53 &
            study.week <= 105 ~ 1,
          
          # post-treatment 2
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week > 105 ~ 0,
          
          # 1 and 4
          # pre-treatment
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week <= 53 ~ 0,
          
          # post-treatment 1
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week >= 54 &
            study.week <= 106 ~ 1,
          
          # post-treatment 2
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week > 106 ~ 0
          
        ),
        
        post2 = case_when(
          
          # 2 and 3
          # pre-treatment
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week < 53 ~ 0,
          
          # post-treatment 1
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week >= 53 &
            study.week <= 105 ~ 0,
          
          # post-treatment 2
          Site %in% c("2A", "2B", "2C", "3A", "3B", "3C") &
            study.week > 105 ~ 1,
          
          # 1 and 4
          # pre-treatment
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week <= 53 ~ 0,
          
          # post-treatment 1
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week >= 54 &
            study.week <= 106 ~ 0,
          
          # post-treatment 2
          Site %in% c("1A", "1B", "1C", "4A", "4B", "4C") &
            study.week > 106 ~ 1
          
        ),
        
        # assign treatments
        ret = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ 0,
                        Site %in% c("1A", "2B", "3B", "4A") ~ 1,
                        Site %in% c("1B", "2A", "3A", "4B") ~ 0),
        
        pil = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ 0,
                        Site %in% c("1A", "2B", "3B", "4A") ~ 0,
                        Site %in% c("1B", "2A", "3A", "4B") ~ 1)
        
      )
    
  )
  
  # return
  return(x.df)
  
}

# 1. calculate complete hazard
# 2. simulate draws
# 3. truncate at the first > 0
sim_lifetime <- function (x) {
  
  # x is each individual deployment
  
  # 1. calculate complete hazard
  calc_spline_pred <- function (x, y) {
    
    # y is the focal iteration
    # extract spline weights w
    w <- t(as.matrix(as.numeric(c(y[paste0("w", "[", x$sex_forest[1], ", 1]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 2]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 3]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 4]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 5]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 6]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 7]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 8]")],
                                  y[paste0("w", "[", x$sex_forest[1], ", 9]")]))))
    
    # multiply and sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
    
    # add to intercept and exponentiate
    blh = as.numeric(exp(y[paste0("a0[", x$sex_forest, "]")] + w.by.b.sum[x$week]))
    
    # total hazard ratio prediction
    hr = exp(log(y$hr_bci) * x$BCI.s +
               log(y$hr_bci_study_week) * x$study.week.s * x$BCI.s +
               log(y$hr_dm) * x$p.dm.s +
               log(y$hr_open) * x$p.open.s +
               log(y$hr_pil_total1) * x$post1 * x$pil +
               log(y$hr_pil_total2) * x$post2 * x$pil +
               log(y$hr_ret_total1) * x$post1 * x$ret +
               log(y$hr_ret_total2) * x$post2 * x$ret)
    
    # full hazard
    full.haz = blh * hr
    
    # return
    return(full.haz)
    
  }
  
  # add to x
  x$haz.1 <- calc_spline_pred(x, iter.draw1)
  x$haz.2 <- calc_spline_pred(x, iter.draw2)
  x$haz.3 <- calc_spline_pred(x, iter.draw3)
  
  # 2. simulate draws
  draw_pois <- function (x) {
    
    x.1 <- x %>%
      
      mutate(
        
        draw.1 = rpois(n(), haz.1),
        draw.2 = rpois(n(), haz.2),
        draw.3 = rpois(n(), haz.3)
        
      )
    
    return(x.1)
    
  }
  
  x.draws <- draw_pois(x)
  
  # 3. truncate at first event, extract lifetime
  x.2 <- data.frame(
    
    deploy = x.draws$deploy[1],
    
    lifetime1 = ifelse(
      
      sum(x.draws$draw.1) > 0,
      which(x.draws$draw.1 > 0)[1],
      max(x.draws$t)
      
    ),
    
    lifetime2 = ifelse(
      
      sum(x.draws$draw.2) > 0,
      which(x.draws$draw.2 > 0)[1],
      max(x.draws$t)
      
    ),
    
    lifetime3 = ifelse(
      
      sum(x.draws$draw.3) > 0,
      which(x.draws$draw.3 > 0)[1],
      max(x.draws$t)
      
    )
    
  )
  
  return(x.2)
  
}

# calculate Kaplan-Meier curve from lifetimes
kap_mei_lifetime <- function (x) {
  
  # x is the data.frame of lifetimes
  
  # function that takes each scenario's lifetimes and determines who's alive at each point
  lifetime_to_survivorship <- function (y) {
    
    prop.alive <- vector(length = 150)
    
    for (i in 1:150) {
      
      prop.alive[i] <- sum(y >= i) / length(y)
      
    }
    
    return(prop.alive)
    
  }
  
  # data.frame to hold it all
  survivorship.all <- data.frame(
    
    t = rep(1:150, times = 3),
    
    S = c(lifetime_to_survivorship(x$lifetime1),
          lifetime_to_survivorship(x$lifetime2),
          lifetime_to_survivorship(x$lifetime3)),
    
    scen = rep(c(1:3), each = 150)
    
  )
  
  # return
  return(survivorship.all)
  
}

#_______________________________________________________________________________________________
# 6b. Loop ----

# split deployments
deploy.split <- split(fates.2, fates.2$deployment.1)

#_______________________________________________________________________________________________

# loop through draws
S.df.all <- data.frame()

for (i in 1:nrow(model.fit.1)) {
  
  # extract focal iteration
  iter.draw1 <- model.fit.1[i, ]
  iter.draw2 <- model.fit.2[i, ]
  iter.draw3 <- model.fit.3[i, ]
  
  # 1. expand deployment dfs - expand_time
  expanded.deploys <- lapply(deploy.split, expand_time, max.t = 150)
  
  # 2. calculate hazard at each week
  # 3. take Poisson draws for each week
  # 4. truncate at the first > 0
  lifetimes <- lapply(expanded.deploys, sim_lifetime)
  
  # bind lifetimes together
  lifetimes.df <- do.call(rbind, lifetimes)
  
  # calculate survivorship
  S.df <- kap_mei_lifetime(lifetimes.df)
  
  # add iter column
  S.df$iter <- i
  
  # bind in (forgot this!)
  S.df.all <- rbind(S.df.all, S.df)
  
  # print status and save to .csv every 50 iterations
  if (i %% 50 == 0) {
    
    print(paste0("Completed iteration ", i, " of 3000"))
    
    write.csv(S.df.all, "S_all_i")
    
  }
  
}

# ~ 13 s per iteration
# 10 h
(13 * 3000) / 3600

#_______________________________________________________________________________________________
# 7. Read data back in ----
#_______________________________________________________________________________________________

S.df <- read.csv("PPCs/S_all_i.csv")

#_______________________________________________________________________________________________
# 8. Final survivorship curves for PPC ----
#_______________________________________________________________________________________________

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

km.ppc <- km_pred(S.df)  

#_______________________________________________________________________________________________
# 8c. Plot with empirical curve ----

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
  
  # simulated prediction envelopes
  geom_line(data = km.ppc,
            aes(x = t,
                y = med),
            color = "aquamarine4",
            linewidth = 1) +
  
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
              fill = NA,
              linetype = "dashed") +
  
  # empirical curve
  geom_line(data = km.df.all,
            aes(x = t,
                y = S),
            linetype = "dashed") +
  
  xlab("Weeks") +
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(4.25, 70))

# 900 x 343
