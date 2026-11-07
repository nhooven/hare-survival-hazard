# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04 - Prepare data for modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 25 Feb 2026
# Date completed: 25 Feb 2026
# Date last modified: 04 Nov 2026
# R version: 4.5.2

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)            # construct basis functions for cyclic splines

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

fates <- read.csv("data_cleaned/fates_forModel.csv")
fates.deploy <- readRDS("data_cleaned/fates_deploy.rds")
day.lookup <- read.csv("data_cleaned/day_lookup.csv")

#_______________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----

# we'll use smooth2random here to make sure everything is parameterized right in the model

#_______________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

# construct raw smooth
smooth.raw <- smoothCon(s(week,
                          bs = "cp",        # cyclic penalized - what I really want
                          k = n.knots),
                        data = fates,
                        absorb.cons = T)    # identifiability constraint

# smooth2random
smooth.rand <- smooth2random(smooth.raw[[1]], 
                             vnames = "",
                             type = 2) # fixed and random

# random effect basis design
re.basis <- smooth.rand$rand$Xr

#_______________________________________________________________________________
# 4. Examine correlations between covariates ----
#_______________________________________________________________________________

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

#_______________________________________________________________________________
# 5. Constant list ----
#_______________________________________________________________________________

constant.list = list(
  
  # dimensions
  N = nrow(fates),              
  nclust = 4, 
  nbasis = ncol(re.basis),
  
  # indices
  cluster = fates$cluster,
  sex = fates$sex + 1,
  deployment = fates$deployment.1,
  
  # spline
  basis = re.basis,
  
  # BCI imputation
  bci.mean = fates.deploy$bci.mean,
  bci.sd = fates.deploy$bci.sd
  
)

#_______________________________________________________________________________
# 6. Data list ----
#_______________________________________________________________________________

# data
data.list <- list(
  
  # responses, split by scenario
  y = fates$y.mort1,
  
  # individual-week covariates
  # treatment
  post1 = fates$post1,
  post2 = fates$post2,
  ret = fates$ret,
  pil = fates$pil,
  
  # collar
  collar = fates$collar,
  
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

#_______________________________________________________________________________
# 7. Save lists ----
#_______________________________________________________________________________

saveRDS(constant.list, "data_cleaned/constants.rds")
saveRDS(data.list, "data_cleaned/data.rds")
