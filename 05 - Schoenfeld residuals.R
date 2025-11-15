# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 05 - Schoenfeld residuals
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 
# Date last modified: 14 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- read.csv("Model outputs/model_1.csv")

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
#_______________________________________________________________________________________________

model.fit.1.hr <- model.fit.1 %>%
  
  dplyr::select(
    
    hr_bci,
    hr_pil_total1,
    hr_pil_total2,
    hr_ret_total1,
    hr_ret_total2
    
  )

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
        study.week,
        study.year.week
        
        ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1))

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
  response,     # which variable should we look for 1s in?
  ci = 0.90
  
  ) {
  
  # only proceed if there are > 0 events
  if (1 %in% x[[response]]) {
    
    # calculate mean risk-weighted covariate value for each posterior HR draw
    # create blank vectors
    rw.cov.vec.mean <- vector(length = nrow(model.fit.1.hr))    # mean
    rw.cov.vec.var <- vector(length = nrow(model.fit.1.hr))    # variance
    
    # loop through all iterations
    for (i in 1:nrow(model.fit.1.hr)) {
  
        rw.cov.vec.mean[i] <- mean(x$BCI.s * model.fit.1.hr$hr_bci[i]) 
        rw.cov.vec.var[i] <- var(x$BCI.s * model.fit.1.hr$hr_bci[i]) 
      
    }
    
    # calculate scaled Schoenfeld residual for all events at risk-time (for each draw)
    # subset of events
    x.event <- x[x[[response]] == 1, ]
    
    # blank matrix
    schoen.matrix <- matrix(data = NA,
                            nrow = nrow(model.fit.1.hr),
                            ncol = sum(x[[response]]))
    
    for (j in 1:sum(x[[response]])) {
      
      schoen.matrix[ ,j] <- (x.event$BCI.s[j] - rw.cov.vec) * (1 / rw.cov.vec.var)
      
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
    
  }
  
  return(schoen.df)
  
}

# apply function
schoen.bci.test <- do.call(rbind, 
                           lapply(split(fates.1, fates.1$study.week), 
                                  schoen_resid_bci,
                                  response = "y.mort.scen1",
                                  ci = 0.90))

# plot test
ggplot(data = schoen.bci.test) +
  
  theme_bw() +
  
  geom_point(aes(x = study.week,
                  y = med)) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_smooth(aes(x = study.week,
                  y = med),
              method = "loess",
              se = T)

# looks like we have a general trend over time. Interacting BCI with study.week could de-trend this!
