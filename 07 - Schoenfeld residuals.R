# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 07 - Schoenfeld residuals
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 17 Nov 2025 
# Date last modified: 01 May 2026 
# R version: 4.4.3

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

# model samples
model.fit.1 <- readRDS("models/model_1.rds")
model.fit.2 <- readRDS("models/model_2.rds")
model.fit.3 <- readRDS("models/model_3.rds")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# BCI by deployment
fates.deploy <- readRDS("Cleaned data/fates_deploy.rds")
fates.deploy.df <- fates.deploy$fates.deploy

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________
# 3a. Betas ----

# function to turn mcmc.list into df and select only columns we need
extract_b <- function (x) {
  
  x.1 <- do.call(rbind, x)
  
  x.1.hr <- as.data.frame(x.1) %>%
    
    dplyr::select(
      
      b_bci,
      b_dm,
      b_open,
      b_post1,
      b_post2,
      b_ret_post1,
      b_ret_post2,
      b_pil_post1,
      b_pil_post2
      
    )
  
  return(x.1.hr)
  
}

#_______________________________________________________________________________

model.fit.1.b <- extract_b(model.fit.1)
model.fit.2.b <- extract_b(model.fit.2)
model.fit.3.b <- extract_b(model.fit.3)

#_______________________________________________________________________________
# 3b. Attribute study-week to fates df ----
#_______________________________________________________________________________

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
  
  # join in BCI data
  left_join(fates.deploy.df %>% dplyr::select(deployment.1, BCI)) %>%
  
  # impute mean BCI 
  replace_na(list(BCI = fates.deploy$bci.mean)) %>%
  
  # standardize variables
  mutate(BCI.s = (BCI - fates.deploy$bci.mean) / fates.deploy$bci.sd,
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________
# 3c. Split dataset by study year ----
#_______________________________________________________________________________

fates.1.pre <- fates.1 %>% filter(post1 == 0 & post2 == 0)
fates.1.post1 <- fates.1 %>% filter(post1 == 1 & post2 == 0)
fates.1.post2 <- fates.1 %>% filter(post1 == 0 & post2 == 1)

#_______________________________________________________________________________
# 4. Write functions to calculate Schoenfeld residuals ----
#_______________________________________________________________________________
# 4a. BCI ----

# we'll calculate this over the entire dataset, because the model is naive to this changing
# over time!

#_______________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_bci <- function (
  
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  ci = 0.90,
  scenario = 1
  
  ) {
  
  # which model fit to use?
  model.fit.b <- case_when(scenario == 1 ~ model.fit.1.b,
                           scenario == 2 ~ model.fit.2.b,
                           scenario == 3 ~ model.fit.3.b)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior b draw
    # create blank vectors
    rw.cov.vec.mean <- vector(length = nrow(model.fit.b))    # mean
    rw.cov.vec.var <- vector(length = nrow(model.fit.b))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.b)) {
  
        rw.cov.vec.mean[i] <- mean(x$BCI.s * model.fit.b$b_bci[i])
        rw.cov.vec.var[i] <- var(x$BCI.s * model.fit.b$b_bci[i])
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrix
    schoen.matrix <- matrix(data = NA,
                            nrow = nrow(model.fit.b),
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
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix[ ,k], ci = ci))[3],
        sd = sd(schoen.matrix[ , k])
        
      )
      
      # bind in
      schoen.df <- rbind(schoen.df, focal.df)
      
    }
    
    return(schoen.df)
    
  }
  
}

#_______________________________________________________________________________
# 4b. Treatment ----
#_______________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_trt <- function (
    
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  post = 1,                      # which period?
  ci = 0.90,
  scenario = 1
  
) {
  
  # which model fit to use?
  model.fit.b <- case_when(scenario == 1 ~ model.fit.1.b,
                           scenario == 2 ~ model.fit.2.b,
                           scenario == 3 ~ model.fit.3.b)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean.ret <- vector(length = nrow(model.fit.b))    # mean
    rw.cov.vec.mean.pil <- vector(length = nrow(model.fit.b))    # mean
    rw.cov.vec.var.ret <- vector(length = nrow(model.fit.b))    # variance
    rw.cov.vec.var.pil <- vector(length = nrow(model.fit.b))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.1.b)) {
      
      # which period?
      
      if (post == 1) {
        
        rw.cov.vec.mean.ret[i] <- mean(model.fit.b$b_post1[i] + x$ret * model.fit.b$b_ret_post1[i]) 
        rw.cov.vec.mean.pil[i] <- mean(model.fit.b$b_post1[i] + x$pil * model.fit.b$b_pil_post1[i]) 
        rw.cov.vec.var.ret[i] <- var(model.fit.b$b_post1[i] + x$ret * model.fit.b$b_ret_post1[i]) 
        rw.cov.vec.var.pil[i] <- var(model.fit.b$b_post1[i] + x$pil * model.fit.b$b_pil_post1[i]) 
        
      }
      
      if (post == 2) {
        
        rw.cov.vec.mean.ret[i] <- mean(model.fit.b$b_post2[i] + x$ret * model.fit.b$b_ret_post2[i]) 
        rw.cov.vec.mean.pil[i] <- mean(model.fit.b$b_post2[i] + x$pil * model.fit.b$b_pil_post2[i]) 
        rw.cov.vec.var.ret[i] <- var(model.fit.b$b_post2[i] + x$ret * model.fit.b$b_ret_post2[i]) 
        rw.cov.vec.var.pil[i] <- var(model.fit.b$b_post2[i] + x$pil * model.fit.b$b_pil_post2[i]) 
        
      }
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrices
    schoen.matrix.ret <- matrix(data = NA,
                                nrow = nrow(model.fit.b),
                                ncol = sum(x[[response]]))
    
    schoen.matrix.pil <- matrix(data = NA,
                                nrow = nrow(model.fit.b),
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
        sd = sd(schoen.matrix.ret[ ,k]),
        post = post,
        trt = "ret"
        
      )
      
      focal.df.pil <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.pil[ ,k], ci = ci))[2],
        med = median(schoen.matrix.pil[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.pil[ ,k], ci = ci))[3],
        sd = sd(schoen.matrix.pil[ ,k]),
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

#_______________________________________________________________________________
# 4c. Landscape variables ----
#_______________________________________________________________________________

# input is a df split into a list by study.week

schoen_resid_lsm <- function (
    
  x,
  response = "y.mort.scen1",     # which variable should we look for 1s in?
  ci = 0.90,
  scenario = 1
  
) {
  
  # which model fit to use?
  model.fit.b <- case_when(scenario == 1 ~ model.fit.1.b,
                           scenario == 2 ~ model.fit.2.b,
                           scenario == 3 ~ model.fit.3.b)
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean.dm <- vector(length = nrow(model.fit.b))    # mean
    rw.cov.vec.mean.open <- vector(length = nrow(model.fit.b))    # mean
    rw.cov.vec.var.dm <- vector(length = nrow(model.fit.b))    # variance
    rw.cov.vec.var.open <- vector(length = nrow(model.fit.b))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.1.b)) {
        
      rw.cov.vec.mean.dm[i] <- mean(x$p.dm.s * model.fit.b$b_dm[i]) 
      rw.cov.vec.mean.open[i] <- mean(x$p.open.s * model.fit.b$b_open[i]) 
      rw.cov.vec.var.dm[i] <- var(x$p.dm.s * model.fit.b$b_dm[i]) 
      rw.cov.vec.var.open[i] <- var(x$p.open.s * model.fit.b$b_open[i]) 
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrices
    schoen.matrix.dm <- matrix(data = NA,
                                nrow = nrow(model.fit.b),
                                ncol = sum(x[[response]]))
    
    schoen.matrix.open <- matrix(data = NA,
                                 nrow = nrow(model.fit.b),
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
        sd = sd(schoen.matrix.dm[ ,k]),
        u.ci = as.numeric(bayestestR::hdi(schoen.matrix.dm[ ,k], ci = ci))[3],
        ls = "dm"
        
      )
      
      focal.df.open <- data.frame(
        
        study.week = x$study.week[1],
        l.ci = as.numeric(bayestestR::hdi(schoen.matrix.open[ ,k], ci = ci))[2],
        med = median(schoen.matrix.open[ ,k]),
        sd = sd(schoen.matrix.open[ ,k]),
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

#_______________________________________________________________________________
# 5. Calculate Schoenfeld residals ----
#_______________________________________________________________________________
# 5a. BCI ----
#_______________________________________________________________________________

# model 1
schoen.bci.1 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_bci,
                               response = "y.mort.scen1",
                               ci = 0.95,
                               scenario = 1))

# model 2
schoen.bci.2 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_bci,
                               response = "y.mort.scen2",
                               ci = 0.95,
                               scenario = 2))

# model 3
schoen.bci.3 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_bci,
                               response = "y.mort.scen3",
                               ci = 0.95,
                               scenario = 3))

#_______________________________________________________________________________
# 5b. Treatment ----
#_______________________________________________________________________________

# model 1
schoen.trt.all.1 <- rbind(
  
  do.call(rbind, 
          lapply(split(fates.1.post1, fates.1.post1$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen1",
                 post = 1,
                 ci = 0.95,
                 scenario = 1)),
  do.call(rbind, 
          lapply(split(fates.1.post2, fates.1.post2$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen1",
                 post = 2,
                 ci = 0.95,
                 scenario = 1))
  
)

# model 2
schoen.trt.all.2 <- rbind(
  
  do.call(rbind, 
          lapply(split(fates.1.post1, fates.1.post1$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen2",
                 post = 1,
                 ci = 0.95,
                 scenario = 2)),
  do.call(rbind, 
          lapply(split(fates.1.post2, fates.1.post2$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen2",
                 post = 2,
                 ci = 0.95,
                 scenario = 2))
  
)

# model 3
schoen.trt.all.3 <- rbind(
  
  do.call(rbind, 
          lapply(split(fates.1.post1, fates.1.post1$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen3",
                 post = 1,
                 ci = 0.95,
                 scenario = 3)),
  do.call(rbind, 
          lapply(split(fates.1.post2, fates.1.post2$study.week), 
                 schoen_resid_trt,
                 response = "y.mort.scen3",
                 post = 2,
                 ci = 0.95,
                 scenario = 3))
  
)

#_______________________________________________________________________________
# 5c. Landscape ----
#_______________________________________________________________________________

# model 1
schoen.lsm.1 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_lsm,
                               response = "y.mort.scen1",
                               ci = 0.95,
                               scenario = 1))

# model 2
schoen.lsm.2 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_lsm,
                               response = "y.mort.scen2",
                               ci = 0.95,
                               scenario = 2))

# model 3
schoen.lsm.3 <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               schoen_resid_lsm,
                               response = "y.mort.scen3",
                               ci = 0.95,
                               scenario = 3))

#_______________________________________________________________________________
# 6. Plot Schoenfeld residals ----
#_______________________________________________________________________________
# 6a. BCI ----
#_______________________________________________________________________________

# bind together
schoen.bci.all <- rbind(schoen.bci.1 %>% mutate(scen = "Scenario 1"),
                        schoen.bci.2 %>% mutate(scen = "Scenario 2"),
                        schoen.bci.3 %>% mutate(scen = "Scenario 3"))


# plot
ggplot(data = schoen.bci.all) +
  
  theme_bw() +
  
  facet_wrap(~ scen) +
  
  geom_point(aes(x = study.week,
                 y = med),
             shape = 21,
             size = 0.7,
             fill=  "white") +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  
  xlab("Study week") +
  ylab("Scaled Schoenfeld residual")

# 667 x 291

#_______________________________________________________________________________
# 6b. Treatment ----
#_______________________________________________________________________________

# bind together
schoen.trt.all <- rbind(schoen.trt.all.1 %>% mutate(scen = "Scenario 1"),
                        schoen.trt.all.2 %>% mutate(scen = "Scenario 2"),
                        schoen.trt.all.3 %>% mutate(scen = "Scenario 3"))


# plot
ggplot(data = schoen.trt.all) +
  
  theme_bw() +
  
  facet_grid(trt ~ scen) +
  
  geom_point(aes(x = study.week,
                 y = med),
             shape = 21,
             size = 0.7,
             fill=  "white") +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  # POST1
  geom_smooth(data = schoen.trt.all %>% filter(post == 1),
              aes(x = study.week,
                  y = med,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  # POST2
  geom_smooth(data = schoen.trt.all %>% filter(post == 2),
              aes(x = study.week,
                  y = med,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  
  xlab("Study week") +
  ylab("Scaled Schoenfeld residual")

# 667 x 470

#_______________________________________________________________________________
# 6c. Landscape ----
#_______________________________________________________________________________

# bind together
schoen.lsm.all <- rbind(schoen.lsm.1 %>% mutate(scen = "Scenario 1"),
                        schoen.lsm.2 %>% mutate(scen = "Scenario 2"),
                        schoen.lsm.3 %>% mutate(scen = "Scenario 3"))


# plot
ggplot(data = schoen.lsm.all) +
  
  theme_bw() +
  
  facet_grid(ls ~ scen) +
  
  geom_point(aes(x = study.week,
                 y = med),
             shape = 21,
             size = 0.7,
             fill=  "white") +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  
  xlab("Study week") +
  ylab("Scaled Schoenfeld residual")

# 667 x 470
