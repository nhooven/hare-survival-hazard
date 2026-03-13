# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02 - Body condition index
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 13 Nov 2025
# Date completed: 13 Nov 2025
# Date last modified: 09 Mar 2026 
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
# 6. Set up deployment-specific dataset for model imputation ----
#_______________________________________________________________________________________________

fates.deploy <- fates.1 %>%
  
  dplyr::select(deployment.1,
                Sex.1,
                HFL,
                Final.mass) %>%
  
  # calculate BCI (implicitly for complete cases)
  mutate(BCI = (Final.mass * 1000) / (HFL * 10))

#_______________________________________________________________________________________________
# 7. Write to file ----
#_______________________________________________________________________________________________

fates.deploy.list <- list(
  
  # deployment df
  fates.deploy = fates.deploy,
  
  # standardization parameters
  bci.mean = mean(fates.deploy$BCI, na.rm = T),
  bci.sd = sd(fates.deploy$BCI, na.rm = T)
  
)

saveRDS(fates.deploy.list, "Cleaned data/fates_deploy.rds")
