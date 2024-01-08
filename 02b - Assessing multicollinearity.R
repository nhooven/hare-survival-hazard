# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 02b - Assessing multicollinearity
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 08 Jan 2024
# Date completed: 
# Date last modified: 08 Jan 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("fates_cleaned.csv")

# keep only columns we need for modeling
fates.1 <- fates %>% dplyr::select(cluster,
                                   Site,
                                   Ear.tag,
                                   Collar.type,
                                   week,
                                   status.num,
                                   Sex.1,
                                   Mass.1,
                                   HFL.1,
                                   Treatment.Retention,
                                   Treatment.Piling)

#_______________________________________________________________________________________________
# 3. Examine potential multicollinearity ----
#_______________________________________________________________________________________________
# 3a. Relationship between sex and morphometry ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # points per sex
  geom_point(data = fates.1,
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
  geom_smooth(data = fates.1 %>% filter(Sex.1 == 0),
              aes(x = HFL.1,
                  y = Mass.1),
              method = "lm",
              color = "red",
              fill = "red",
              alpha = 0.25) +
  
  # male
  geom_smooth(data = fates.1 %>% filter(Sex.1 == 1),
              aes(x = HFL.1,
                  y = Mass.1),
              method = "lm",
              color = "darkgray",
              fill = "darkgray",
              alpha = 0.25)

# here the female intercept and slope is higher

# let's run the regression out of curiosity
sex.reg <- lm(data = fates.1,
              Mass.1 ~ Sex.1 + HFL.1)

summary(sex.reg)



# correlation between continuous covariates (mass and HFL)
cor(x = fates.1[ ,c("Mass.1", "HFL.1")], method = "pearson")

ggplot(data = fates.1,
       aes(x = Mass.1,
           y = HFL.1)) +
  theme_bw() +
  geom_point()

# kruskal-wallis tests for continuous by categorical
# mass by sex
kruskal.test(x = fates.1$Mass.1, g = fates.1$Sex.1)

ggplot(data = fates.1,
       aes(fill = as.factor(Sex.1),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by sex
kruskal.test(x = fates.1$HFL.1, g = fates.1$Sex.1)

ggplot(data = fates.1,
       aes(fill = as.factor(Sex.1),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# mass by retention
kruskal.test(x = fates.1$Mass.1, g = fates.1$Treatment.Retention)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Retention),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by retention
kruskal.test(x = fates.1$HFL.1, g = fates.1$Treatment.Retention)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Retention),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# mass by piling
kruskal.test(x = fates.1$Mass.1, g = fates.1$Treatment.Piling)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Piling),
           x = Mass.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# hfl by Piling
kruskal.test(x = fates.1$HFL.1, g = fates.1$Treatment.Piling)

ggplot(data = fates.1,
       aes(fill = as.factor(Treatment.Piling),
           x = HFL.1)) +
  theme_bw() +
  geom_density(alpha = 0.5) +
  theme(legend.position = c(0.8, 0.8))

# variance inflation factors (general Poisson regression)
vif.model <- glm(status.num ~ Sex.1 + Treatment.Retention + Treatment.Piling + Mass.1 + HFL.1,
                 data = fates.1,
                 family = poisson)

car::vif(vif.model)