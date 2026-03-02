# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02 - Body condition index
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Nov 2025
# Date completed: 13 Nov 2025
# Date last modified: 03 Feb 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_final_cleaned.csv")

#_______________________________________________________________________________________________
# 3. Data cleaning ----
#_______________________________________________________________________________________________

# pull out one obs per deployment
fates.1 <- fates %>%
  
  group_by(deployment.1) %>%
  
  slice(1)

# number of missing observations
sum(is.na(fates.1$HFL))
sum(is.na(fates.1$Final.mass))

# keep only complete cases for plotting and mean extraction
fates.cc <- fates.1 %>%
  
  filter(is.na(HFL) == F &
         is.na(Final.mass) == F) %>%
  
  # calculate body condition index (g / mm)
  mutate(BCI = (Final.mass * 1000) / (HFL * 10))

#_______________________________________________________________________________________________
# 4. Complete case plots ----
#_______________________________________________________________________________________________
# 4a. Mass vs HFL relationship ----
#_______________________________________________________________________________________________

ggplot(data = fates.cc,
       aes(x = HFL * 10,
           y = Final.mass * 1000,
           color = as.factor(Sex.1),
           fill = as.factor(Sex.1),
           group = as.factor(Sex.1),
           shape = as.factor(Sex.1))) +
  
  theme_bw() + 
  
  geom_point() +
  
  geom_smooth(method = "lm",
              aes(linetype = as.factor(Sex.1))) +
  
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.75)) +
  
  xlab("Right hind foot length (mm)") +
  ylab("Body mass (g)")

#_______________________________________________________________________________________________
# 4b. Sex-specific mass distributions ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_density(data = fates.cc,
               aes(x = Final.mass * 1000,
                   color = as.factor(Sex.1),
                   fill = as.factor(Sex.1)),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.75)) +
  
  xlab("Body mass (g)")

#_______________________________________________________________________________________________
# 4c. Sex-specific HFL distributions ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_histogram(data = fates.cc,
               aes(x = HFL * 10,
                   color = as.factor(Sex.1),
                   fill = as.factor(Sex.1)),
               alpha = 0.25,
               bins = length(unique(fates.cc$HFL))) +
  
  facet_wrap(~ as.factor(Sex.1)) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  xlab("Right hind foot length (mm)")

#_______________________________________________________________________________________________
# 4d. Sex-specific BCI distributions ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_density(data = fates.cc,
                 aes(x = BCI,
                     color = as.factor(Sex.1),
                     fill = as.factor(Sex.1)),
                 alpha = 0.25) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  scale_x_continuous(breaks = c(7, 8, 9, 10, 11, 12)) +
  
  xlab("Body condition index")

#_______________________________________________________________________________________________
# 5. Calculate sex-specific means for each variable ----
#_______________________________________________________________________________________________

fates.cc.summary <- fates.cc %>%
  
  group_by(Sex.1) %>%
  
  summarize(
    
    mass.mean = mean(Final.mass),
    mass.cv = sd(Final.mass) / mean(Final.mass),
    hfl.mean = mean(HFL),
    hfl.cv = sd(HFL) / mean(HFL),
    bci.mean = mean(BCI),
    bci.cv = sd(BCI) / mean(BCI)
    
    )

fates.cc.summary

#_______________________________________________________________________________________________
# 6. Use each imputation strategy to add BCI values to each missing observation ----

# STRATEGY 1: impute missing mass and HFL independently, then calculate BCI
# STRATEGY 2: impute mean BCI for any incomplete observation (i.e., missing either variable)

#_______________________________________________________________________________________________
# 6a. Write function ----
#_______________________________________________________________________________________________

impute_bci <- function (x) {
  
  # deployment has both an HFL and a mass
  if (is.na(x$HFL[1]) == F & is.na(x$Final.mass[1]) == F) {
    
    # S1 and S2 are the same in this case
    x$BCI.1 = (x$Final.mass * 1000) / (x$HFL * 10)
    x$BCI.2 = (x$Final.mass * 1000) / (x$HFL * 10)
    
  }
  
  # deployment is missing HFL but not mass
  if (is.na(x$HFL[1]) == T & is.na(x$Final.mass[1]) == F) {
    
    # impute missing morphometric
    x$HFL = fates.cc.summary$hfl.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
    # S1
    x$BCI.1 = (x$Final.mass * 1000) / (x$HFL * 10)
    
    # S2
    x$BCI.2 = fates.cc.summary$bci.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
  }
  
  # deployment is missing mass but not HFL
  if (is.na(x$HFL[1]) == F & is.na(x$Final.mass[1]) == T) {
    
    # impute missing morphometric
    x$Final.mass = fates.cc.summary$mass.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
    # S1
    x$BCI.1 = (x$Final.mass * 1000) / (x$HFL * 10)
    
    # S2
    x$BCI.2 = fates.cc.summary$bci.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
  }
  
  # deployment is missing both HFL and mass
  if (is.na(x$HFL[1]) == T & is.na(x$Final.mass[1]) == T) {
    
    # S1
    x$BCI.1 = fates.cc.summary$bci.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
    # S2
    x$BCI.2 = fates.cc.summary$bci.mean[fates.cc.summary$Sex.1 == x$Sex.1[1]]
    
  }
  
  return(x)
  
}

#_______________________________________________________________________________________________
# 6b. Apply function ----
#_______________________________________________________________________________________________

fates.imputed <- do.call(rbind, lapply(split(fates, fates$deployment.1), impute_bci))

#_______________________________________________________________________________________________
# 6c. Compare distributions for each strategy ----
#_______________________________________________________________________________________________

# consolidate
fates.imputed.1 <- fates.imputed %>%
  
  group_by(deployment.1) %>%
  
  slice(1)

# plot a 1:1 graph
ggplot(data = fates.imputed.1,
       aes(x = BCI.1,
           y = BCI.2)) +
  
  theme_bw() +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  geom_point() +
  
  coord_cartesian(xlim = c(6, 13),
                  ylim = c(6, 13))

# sex-specific distributions
ggplot(data = fates.imputed.1 %>% pivot_longer(cols = c("BCI.1", "BCI.2"))) +
  
  theme_bw() +
  
  facet_wrap(~ name) +
  
  geom_density(aes(x = value,
                   color = as.factor(Sex.1),
                   fill = as.factor(Sex.1)),
               alpha = 0.25) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  
  scale_x_continuous(breaks = c(7, 8, 9, 10, 11, 12)) +
  
  xlab("Body condition index")

# unsurprisingly, there is more prob. density at the sex-specific means for BCI.2
# we can try both in the model and see if these imputed data still have a signal

#_______________________________________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________________________________

write.csv(fates.imputed, "Cleaned data/fates_final_cleaned_2.csv")
