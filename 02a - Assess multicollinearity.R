# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02a - Assess multicollinearity
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Jan 2024
# Date completed: 08 Jan 2024
# Date last modified: 28 Sep 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_cleaned_09_28_2024.csv")

#_______________________________________________________________________________________________
# 3. Examine potential multicollinearity ----
#_______________________________________________________________________________________________
# 3a. Relationship between sex and morphometry ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # points per sex
  geom_point(data = fates,
             aes(x = HFL.1,
                 y = Mass.1,
                 color = as.factor(Sex.1),
                 shape = as.factor(Sex.1))) +
  
  scale_color_manual(values = c("red", "darkgray")) +
  
  theme(legend.position = "none") +
  
  # mean intercepts
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  # sex-specific regressions
  # female
  geom_smooth(data = fates %>% filter(Sex.1 == 0),
              aes(x = HFL.1,
                  y = Mass.1),
              method = "lm",
              color = "red",
              fill = "red",
              alpha = 0.25) +
  
  # male
  geom_smooth(data = fates %>% filter(Sex.1 == 1),
              aes(x = HFL.1,
                  y = Mass.1),
              method = "lm",
              color = "darkgray",
              fill = "darkgray",
              alpha = 0.25)

# here the female intercept and slope is higher

# let's run the regression out of curiosity
sex.reg <- lm(data = fates,
              Mass.1 ~ Sex.1 + HFL.1)

summary(sex.reg)

# this is pretty convincing evidence that sex and HFL can reliably predict 
# trends in mass

# correlation between covariates
# mass and HFL
cor.test(x = fates$Mass.1, 
         y = fates$HFL.1,
         method = "pearson")

# this is a somewhat strong correlation

# mass and sex
cor.test(x = fates$Mass.1, 
         y = fates$Sex.1,
         method = "pearson")

# slightly weaker

# hfl and sex
cor.test(x = fates$HFL.1, 
         y = fates$Sex.1,
         method = "pearson")

# very weak



# variance inflation factors (general Poisson regression)
vif.model <- glm(mort ~ Sex.1 + trt.ret + trt.ret + Mass.1 + HFL.1,
                 data = fates,
                 family = poisson)

car::vif(vif.model)

# we see somewhat elevated VIFs for mass and hfl (in the basic Poisson glm)

#_______________________________________________________________________________________________
# 3b. Principal components for hare size ----
#_______________________________________________________________________________________________

# since both variables are z-score standardized, it is trivial to compute the first PC
cor(x = fates[ ,c("Mass.1", "HFL.1", "Sex.1", "PC1")], 
    method = "pearson")

# here each variable is correlated the same to the first PC
# so it captures much of the variability both variables

#_______________________________________________________________________________________________
# 4. Body condition indices ----

# I looked at the residual BCI (from Murray 2002) and it was highly (r = 0.89)
# correlated to mass

#_______________________________________________________________________________________________

# examine correlations
cor(x = fates[ ,c("Mass.1", "HFL.1", "Sex.1", "PC1", "BCI")], 
    method = "pearson")

# TENTATIVE CONCLUSIONS
# Include the body size synthetic variable (PC1)
# Body condition index (BCI)
# And MAYBE sex as it is less strongly correlated to PC1 than mass
