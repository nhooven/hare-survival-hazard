# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04 - Prepare data for modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Feb 2026
# Date completed: 25 Feb 2026
# Date last modified: 09 Mar 2026
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)            # construct basis functions for cyclic splines

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_forModel.csv")
fates.deploy <- readRDS("Cleaned data/fates_deploy.rds")
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# create the basis matrix (cyclic spline)
basis <- cSplineDes(fates$week,
                    knots = knot.list,
                    ord = 4)                # cubic

#_______________________________________________________________________________________________
# 4. Examine correlations between covariates ----
#_______________________________________________________________________________________________

# continuous
covar.contin <- fates %>%
  
  dplyr::select(study.week,
                p.dm,
                p.o)

cor(covar.contin) 
# |r| all < 0.4

# continuous by binary
# we'll do plots here
covar.both <- fates %>%
  
  dplyr::select(post1,
                post2,
                ret,
                pil) %>%
  
  bind_cols(covar.contin) %>%
  
  pivot_longer(cols = c(study.week,
                        p.dm,
                        p.o))

# post1
ggplot(covar.both) +
  
  theme_bw() +
  
  facet_wrap(~ name,
             scales = "free") +
  
  geom_boxplot(aes(x = as.factor(post1),
                   y = value)) +
  
  geom_point(aes(x = as.factor(post1),
                 y = value),
             alpha = 0.5)

# looks fine

# post2
ggplot(covar.both) +
  
  theme_bw() +
  
  facet_wrap(~ name,
             scales = "free") +
  
  geom_boxplot(aes(x = as.factor(post2),
                   y = value)) +
  
  geom_point(aes(x = as.factor(post2),
                 y = value),
             alpha = 0.5)

# ret
ggplot(covar.both) +
  
  theme_bw() +
  
  facet_wrap(~ name,
             scales = "free") +
  
  geom_boxplot(aes(x = as.factor(ret),
                   y = value)) +
  
  geom_point(aes(x = as.factor(ret),
                 y = value),
             alpha = 0.5)

# weird for the proportions, but there are only a few values these can take
# ret
ggplot(covar.both) +
  
  theme_bw() +
  
  facet_wrap(~ name,
             scales = "free") +
  
  geom_boxplot(aes(x = as.factor(pil),
                   y = value)) +
  
  geom_point(aes(x = as.factor(pil),
                 y = value),
             alpha = 0.5)

# not too bad here

#_______________________________________________________________________________________________
# 5. Constant list ----
#_______________________________________________________________________________________________

constant.list = list(
  
  # dimensions
  N = nrow(fates),              
  nclust = 4, 
  nbasis = ncol(basis),
  
  # indices
  cluster = fates$cluster,
  sex = fates$sex + 1,
  deployment = fates$deployment.1,
  
  # spline
  basis = basis,
  
  # BCI imputation
  bci.mean = fates.deploy$bci.mean,
  bci.sd = fates.deploy$bci.sd
  
)

#_______________________________________________________________________________________________
# 6. Data list ----
#_______________________________________________________________________________________________

# data
data.list <- list(
  
  # responses, split by scenario
  y_mort_1 = fates$y.mort.scen1,
  y_mort_2 = fates$y.mort.scen2,
  y_mort_3 = fates$y.mort.scen3,
  
  # individual-week covariates
  # treatment
  post1 = fates$post1,
  post2 = fates$post2,
  ret = fates$ret,
  pil = fates$pil,
  
  # extrinsic
  study_week = (fates$study.week - mean(fates$study.week)) / sd(fates$study.week),
  
  # landscape
  dm = (fates$p.dm - mean(fates$p.dm)) / sd(fates$p.dm),
  open = (fates$p.o - mean(fates$p.o)) / sd(fates$p.o),
  
  # deployment covariates
  dep.sex = fates.deploy$fates.deploy$Sex.1,
  hfl = fates.deploy$fates.deploy$HFL,
  mass = fates.deploy$fates.deploy$Final.mass,
  bci = fates.deploy$fates.deploy$BCI
  
)

#_______________________________________________________________________________________________
# 7. Save lists ----
#_______________________________________________________________________________________________

saveRDS(constant.list, "Cleaned data/constants.rds")
saveRDS(data.list, "Cleaned data/data.rds")
