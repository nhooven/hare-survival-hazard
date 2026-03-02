# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 07 - Schoenfeld residuals
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 17 Nov 2025 
# Date last modified: 27 Feb 2026 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

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
# 3a. Hazard ratios ----

# function to turn mcmc.list into df and select only columns we need
extract_hr <- function (x) {
  
  x.1 <- do.call(rbind, x)
  
  x.1.hr <- as.data.frame(x.1) %>%
    
    dplyr::select(
      
      hr_bci,
      hr_bci_study_week,
      hr_pil_total1,
      hr_pil_total2,
      hr_ret_total1,
      hr_ret_total2,
      hr_dm,
      hr_open
      
    )
  
  return(x.1.hr)
  
}

#_______________________________________________________________________________________________

model.fit.1.hr <- extract_hr(model.fit.1)
model.fit.2.hr <- extract_hr(model.fit.2)
model.fit.3.hr <- extract_hr(model.fit.3)

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
  
  # standardize variables
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1),
         study.week.s = (study.week - mean(study.week)) / sd(study.week),
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________________________
# 3c. Split dataset by study year ----
#_______________________________________________________________________________________________

fates.1.pre <- fates.1 %>% filter(post1 == 0 & post2 == 0)
fates.1.post1 <- fates.1 %>% filter(post1 == 1 & post2 == 0)
fates.1.post2 <- fates.1 %>% filter(post1 == 0 & post2 == 1)

#_______________________________________________________________________________________________
# 4. Write functions to calculate Schoenfeld residuals ----
#_______________________________________________________________________________________________
# 4a. BCI ----

# we'll calculate this over the entire dataset, because the model is naive to this changing
# over time!

#_______________________________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_bci <- function (
  
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  ci = 0.90,
  scenario = 1
  
  ) {
  
  # which model fit to use?
  model.fit.hr <- case_when(scenario == 1 ~ model.fit.1.hr,
                            scenario == 2 ~ model.fit.2.hr,
                            scenario == 3 ~ model.fit.3.hr)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean <- vector(length = nrow(model.fit.hr))    # mean
    rw.cov.vec.var <- vector(length = nrow(model.fit.hr))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.hr)) {
  
        rw.cov.vec.mean[i] <- mean(x$BCI.s * model.fit.hr$hr_bci[i] +
                                   x$BCI.s * x$study.week.s * model.fit.hr$hr_bci_study_week[i]) 
        rw.cov.vec.var[i] <- var(x$BCI.s * model.fit.hr$hr_bci[i] +
                                 x$BCI.s * x$study.week.s * model.fit.hr$hr_bci_study_week[i]) 
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrix
    schoen.matrix <- matrix(data = NA,
                            nrow = nrow(model.fit.hr),
                            ncol = sum(x[[response]]))
    
    for (j in 1:sum(x[[response]])) {
      
      schoen.matrix[ ,j] <- (x.event$BCI.s[j] - rw.cov.vec.mean) * (1 / rw.cov.vec.var)
      
    }
    
    # extract the median and upper/lower credible intervals for each observation
    # blank df
    schoen.df <- data.frame()
    
    for (k in 1:ncol(schoen.matrix)) {
      
      focal.df <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix[ ,k], ci = ci))[2],
        med = median(schoen.matrix[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix[ ,k], ci = ci))[3]
        
      )
      
      # bind in
      schoen.df <- rbind(schoen.df, focal.df)
      
    }
    
    return(schoen.df)
    
  }
  
}

#_______________________________________________________________________________________________
# 4b. Treatment ----
#_______________________________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_trt <- function (
    
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  post = 1,                      # which period?
  ci = 0.90,
  scenario = 1
  
) {
  
  # which model fit to use?
  model.fit.hr <- case_when(scenario == 1 ~ model.fit.1.hr,
                            scenario == 2 ~ model.fit.2.hr,
                            scenario == 3 ~ model.fit.3.hr)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean.ret <- vector(length = nrow(model.fit.hr))    # mean
    rw.cov.vec.mean.pil <- vector(length = nrow(model.fit.hr))    # mean
    rw.cov.vec.var.ret <- vector(length = nrow(model.fit.hr))    # variance
    rw.cov.vec.var.pil <- vector(length = nrow(model.fit.hr))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.1.hr)) {
      
      # which period?
      
      if (post == 1) {
        
        rw.cov.vec.mean.ret[i] <- mean(x$ret * model.fit.hr$hr_ret_total1[i]) 
        rw.cov.vec.mean.pil[i] <- mean(x$pil * model.fit.hr$hr_pil_total1[i]) 
        rw.cov.vec.var.ret[i] <- var(x$ret * model.fit.hr$hr_ret_total1[i]) 
        rw.cov.vec.var.pil[i] <- var(x$pil * model.fit.hr$hr_pil_total1[i]) 
        
      }
      
      if (post == 2) {
        
        rw.cov.vec.mean.ret[i] <- mean(x$ret * model.fit.hr$hr_ret_total2[i]) 
        rw.cov.vec.mean.pil[i] <- mean(x$pil * model.fit.hr$hr_pil_total2[i]) 
        rw.cov.vec.var.ret[i] <- var(x$ret * model.fit.hr$hr_ret_total2[i]) 
        rw.cov.vec.var.pil[i] <- var(x$pil * model.fit.hr$hr_pil_total2[i])
        
      }
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrices
    schoen.matrix.ret <- matrix(data = NA,
                                nrow = nrow(model.fit.hr),
                                ncol = sum(x[[response]]))
    
    schoen.matrix.pil <- matrix(data = NA,
                                nrow = nrow(model.fit.hr),
                                ncol = sum(x[[response]]))
    
    for (j in 1:sum(x[[response]])) {
      
      schoen.matrix.ret[ ,j] <- (x.event$ret[j] - rw.cov.vec.mean.ret) * (1 / rw.cov.vec.var.ret)
      
      schoen.matrix.pil[ ,j] <- (x.event$pil[j] - rw.cov.vec.mean.pil) * (1 / rw.cov.vec.var.pil)
      
    }
    
    # extract the median and upper/lower credible intervals for each observation
    # blank dfs
    schoen.df.ret <- data.frame()
    schoen.df.pil <- data.frame()
    
    for (k in 1:ncol(schoen.matrix.ret)) {
      
      focal.df.ret <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.ret[ ,k], ci = ci))[2],
        med = median(schoen.matrix.ret[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.ret[ ,k], ci = ci))[3],
        post = post,
        trt = "ret"
        
      )
      
      focal.df.pil <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.pil[ ,k], ci = ci))[2],
        med = median(schoen.matrix.pil[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.pil[ ,k], ci = ci))[3],
        post = post,
        trt = "pil"
        
      )
      
      # bind in
      schoen.df.ret <- rbind(schoen.df.ret, focal.df.ret)
      schoen.df.pil <- rbind(schoen.df.pil, focal.df.pil)
      
    }
    
    # bind together
    schoen.df <- rbind(schoen.df.ret, schoen.df.pil)
    
    return(schoen.df)
    
  }
  
}

#_______________________________________________________________________________________________
# 4c. Landscape variables ----
#_______________________________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_lsm <- function (
    
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  ci = 0.90,
  scenario = 1
  
) {
  
  # which model fit to use?
  model.fit.hr <- case_when(scenario == 1 ~ model.fit.1.hr,
                            scenario == 2 ~ model.fit.2.hr,
                            scenario == 3 ~ model.fit.3.hr)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean.dm <- vector(length = nrow(model.fit.hr))    # mean
    rw.cov.vec.mean.open <- vector(length = nrow(model.fit.hr))    # mean
    rw.cov.vec.var.dm <- vector(length = nrow(model.fit.hr))    # variance
    rw.cov.vec.var.open <- vector(length = nrow(model.fit.hr))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.1.hr)) {
        
      rw.cov.vec.mean.dm[i] <- mean(x$p.dm.s * model.fit.hr$hr_dm[i]) 
      rw.cov.vec.mean.open[i] <- mean(x$p.open.s * model.fit.hr$hr_open[i]) 
      rw.cov.vec.var.dm[i] <- var(x$p.dm.s * model.fit.hr$hr_dm[i]) 
      rw.cov.vec.var.open[i] <- var(x$p.open.s * model.fit.hr$hr_open[i]) 
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrices
    schoen.matrix.dm <- matrix(data = NA,
                                nrow = nrow(model.fit.hr),
                                ncol = sum(x[[response]]))
    
    schoen.matrix.open <- matrix(data = NA,
                                 nrow = nrow(model.fit.hr),
                                 ncol = sum(x[[response]]))
    
    for (j in 1:sum(x[[response]])) {
      
      schoen.matrix.dm[ ,j] <- (x.event$p.dm.s[j] - rw.cov.vec.mean.dm) * (1 / rw.cov.vec.var.dm)
      
      schoen.matrix.open[ ,j] <- (x.event$p.open.s[j] - rw.cov.vec.mean.open) * (1 / rw.cov.vec.var.open)
      
    }
    
    # extract the median and upper/lower credible intervals for each observation
    # blank dfs
    schoen.df.dm <- data.frame()
    schoen.df.open <- data.frame()
    
    for (k in 1:ncol(schoen.matrix.dm)) {
      
      focal.df.dm <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.dm[ ,k], ci = ci))[2],
        med = median(schoen.matrix.dm[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.dm[ ,k], ci = ci))[3],
        ls = "dm"
        
      )
      
      focal.df.open <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.open[ ,k], ci = ci))[2],
        med = median(schoen.matrix.open[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.open[ ,k], ci = ci))[3],
        ls = "open"
        
      )
      
      # bind in
      schoen.df.dm <- rbind(schoen.df.dm, focal.df.dm)
      schoen.df.open <- rbind(schoen.df.open, focal.df.open)
      
    }
    
    # bind together
    schoen.df <- rbind(schoen.df.dm, schoen.df.open)
    
    return(schoen.df)
    
  }
  
}

#_______________________________________________________________________________________________
# 5. Calculate Schoenfeld residals ----
#_______________________________________________________________________________________________
# 5a. BCI ----
#_______________________________________________________________________________________________

# model 1
schoen.bci.test.1 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_bci,
                                    response = "y.mort.scen1",
                                    ci = 0.90,
                                    scenario = 1))

# model 2
schoen.bci.test.2 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_bci,
                                    response = "y.mort.scen2",
                                    ci = 0.90,
                                    scenario = 2))

# model 3
schoen.bci.test.3 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_bci,
                                    response = "y.mort.scen3",
                                    ci = 0.90,
                                    scenario = 3))

# plot test
ggplot(data = schoen.bci.test.3) +
  
  theme_bw() +
  
  geom_point(aes(x = study.week,
                  y = med)) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med),
              method = "loess",
              se = T) +
  
  coord_cartesian(ylim = c(-50, 50))

# still some weirdly large residuals between weeks ~20-50
# I think we had very few mortalities here. The rest of the study seems fine

#_______________________________________________________________________________________________
# 5b. Treatment ----
#_______________________________________________________________________________________________

# model 1
# apply function
schoen.trt.post1 <- do.call(rbind, 
                           lapply(split(fates.1.post1, fates.1.post1$study.week), 
                                  schoen_resid_trt,
                                  response = "y.mort.scen1",
                                  post = 1,
                                  ci = 0.90))

schoen.trt.post2 <- do.call(rbind, 
                            lapply(split(fates.1.post2, fates.1.post2$study.week), 
                                   schoen_resid_trt,
                                   response = "y.mort.scen1",
                                   post = 2,
                                   ci = 0.90))

# bind together for plotting
schoen.trt.all <- rbind(schoen.trt.post1, schoen.trt.post2)

# plot test
ggplot(data = schoen.trt.all) +
  
  theme_bw() +
  
  facet_wrap(~ trt) +
  
  geom_point(aes(x = study.week,
                 y = med,
                 color = as.factor(post))) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med,
                  group = as.factor(post)),
              method = "gam",
              se = T) +
  
  theme(legend.position = "none")

# well, looks like we figured out the PH issue!

# model 2
# apply function
schoen.trt.post1.2 <- do.call(rbind, 
                            lapply(split(fates.1.post1, fates.1.post1$study.week), 
                                   schoen_resid_trt,
                                   response = "y.mort.scen2",
                                   post = 1,
                                   ci = 0.90,
                                   scenario = 2))

schoen.trt.post2.2 <- do.call(rbind, 
                            lapply(split(fates.1.post2, fates.1.post2$study.week), 
                                   schoen_resid_trt,
                                   response = "y.mort.scen2",
                                   post = 2,
                                   ci = 0.90,
                                   scenario = 2))

# bind together for plotting
schoen.trt.all.2 <- rbind(schoen.trt.post1.2, schoen.trt.post2.2)

# plot test
ggplot(data = schoen.trt.all.2) +
  
  theme_bw() +
  
  facet_wrap(~ trt) +
  
  geom_point(aes(x = study.week,
                 y = med,
                 color = as.factor(post))) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med,
                  group = as.factor(post)),
              method = "gam",
              se = T) +
  
  theme(legend.position = "none")

# model 3
# apply function
schoen.trt.post1.3 <- do.call(rbind, 
                              lapply(split(fates.1.post1, fates.1.post1$study.week), 
                                     schoen_resid_trt,
                                     response = "y.mort.scen3",
                                     post = 1,
                                     ci = 0.90,
                                     scenario = 3))

schoen.trt.post2.3 <- do.call(rbind, 
                              lapply(split(fates.1.post2, fates.1.post2$study.week), 
                                     schoen_resid_trt,
                                     response = "y.mort.scen3",
                                     post = 2,
                                     ci = 0.90,
                                     scenario = 3))

# bind together for plotting
schoen.trt.all.3 <- rbind(schoen.trt.post1.3, schoen.trt.post2.3)

# plot test
ggplot(data = schoen.trt.all.3) +
  
  theme_bw() +
  
  facet_wrap(~ trt) +
  
  geom_point(aes(x = study.week,
                 y = med,
                 color = as.factor(post))) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med,
                  group = as.factor(post)),
              method = "gam",
              se = T) +
  
  theme(legend.position = "none")

#_______________________________________________________________________________________________
# 5c. Landscape ----
#_______________________________________________________________________________________________

# model 1
schoen.lsm.test.1 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_lsm,
                                    response = "y.mort.scen1",
                                    ci = 0.90,
                                    scenario = 1))

# model 2
schoen.lsm.test.2 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_lsm,
                                    response = "y.mort.scen2",
                                    ci = 0.90,
                                    scenario = 2))

# model 3
schoen.lsm.test.3 <- do.call(rbind, 
                             lapply(split(fates.1, fates.1$study.week), 
                                    schoen_resid_lsm,
                                    response = "y.mort.scen3",
                                    ci = 0.90,
                                    scenario = 3))

# plot test
ggplot(data = schoen.lsm.test.3) +
  
  theme_bw() +
  
  facet_wrap(~ ls) +
  
  geom_point(aes(x = study.week,
                 y = med)) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med),
              method = "gam",
              se = T)

# exactly what I hoped to see!
