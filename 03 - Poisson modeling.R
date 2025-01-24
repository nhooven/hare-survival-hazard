# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03b - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 29 Dec 2023
# Date last modified: 22 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(mgcv)            # construct basis functions for cyclic splines

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_cleaned_01_02_2025.csv")

# keep only columns we need for modeling
fates.1 <- fates %>% dplyr::select(cluster,
                                   Site,
                                   Ear.tag,
                                   week,
                                   year,
                                   mort.pred,
                                   mort.oth,
                                   cens,
                                   Collar.type.1,
                                   Sex.1,
                                   Mass.1,
                                   HFL.1,
                                   PC1,
                                   BCI.1,
                                   post.trt,
                                   trt.ret,
                                   trt.pil) %>%
  
  # create "mort" variable (since most morts are predation or unknown)
  mutate(mort = ifelse(mort.pred == 1 | mort.oth == 1,
                       1,
                       0))

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 7 + 1

knot.list <- quantile(fates.1$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# create the basis matrix (cyclic spline)
basis <- cSplineDes(fates.1$week,
                    knots = knot.list,
                    ord = 4)                # cubic

#_______________________________________________________________________________________________
# 4. Examine correlation between variables ----
#_______________________________________________________________________________________________
# 4a. Continuous covariates ----
#_______________________________________________________________________________________________

round(cor(fates[ , c("Mass.1", "HFL.1", "PC1", "BCI.1")]), digits = 3)

plot(fates[ , c("Mass.1", "HFL.1", "PC1", "BCI.1")])

# yep, these are strange, definitely issues with orthogonality if we include BCI

# we'll start with mass and HFL

#_______________________________________________________________________________________________
# 4b. Sex and morphometrics ----
#_______________________________________________________________________________________________

kruskal.test(Mass.1 ~ Sex.1, data = fates)    # different mass by sex
mean(fates$Mass.1[fates$Sex.1 == 1])          # male
mean(fates$Mass.1[fates$Sex.1 == 0])          # female

kruskal.test(HFL.1 ~ Sex.1, data = fates)     # not different HFL by sex

plot(fates[ , c("Mass.1", "HFL.1", "Sex.1")])

# these are statistically different, but are they biologically meaningful? 

#_______________________________________________________________________________________________
# 4c. Treatment and morphometrics ----
#_______________________________________________________________________________________________

kruskal.test(Mass.1 ~ post.trt, data = fates)    # slightly different mass by treatment status
mean(fates$Mass.1[fates$post.trt == 0])          # pre
mean(fates$Mass.1[fates$post.trt == 1])          # post

kruskal.test(HFL.1 ~ post.trt, data = fates)     # different HFL by treatment status

plot(fates[ , c("Mass.1", "HFL.1", "post.trt")])

# difficult to test differences with interactions, let's keep them for now

#_______________________________________________________________________________________________
# 5. Build data list ----
#_______________________________________________________________________________________________

fates.stan.1 <- list(N = nrow(fates.1),
                     y_mort = fates.1$mort,
                     y_cens = fates.1$cens,
                     t = fates.1$week,
                     collar = fates.1$Collar.type.1,
                     sex = fates.1$Sex.1,
                     mas = fates.1$Mass.1,
                     hfl = fates.1$HFL.1,
                     trt = fates.1$post.trt,
                     ret = fates.1$trt.ret,
                     pil = fates.1$trt.pil,
                     clust = fates.1$cluster,
                     nclust = 4,
                     basis = basis,
                     nbasis = ncol(basis))

#_______________________________________________________________________________________________
# 6. Fit model ----
#_______________________________________________________________________________________________

hazard_model <- rstan::stan(
  
  file = "hazard_model.stan",
  data = fates.stan.1,
  chains = 1,
  warmup = 1000,
  iter = 2000
  
)

print(hazard_model, pars = c("b_ret",
                             "b_trt_ret",
                             "b_pil",
                             "b_trt_pil"))

plot(hazard_model, pars = c("b_trt_ret"))

traceplot(hazard_model, pars = c("a0_mean"))

# 01-23-2025: Not really any evidence that the censoring rate differs between clusters
# let's full pool this

# penalized spline coefficients look pretty good, definitely some variability between sites

# and now - for the coefficients!

# coefs look good! let's calculate the total treatments
# ret
data.frame(wrabbit = round(exp(c(-0.54, -0.54 + 0.58)), digits = 3),
           crazy = round(exp(c(-0.63, -0.63 + 0.48)), digits = 3),
           juice = round(exp(c(-0.70, -0.70 + 0.47)), digits = 3),
           chop = round(exp(c(-0.70, -0.70 + 0.40)), digits = 3))

# ret
data.frame(wrabbit = round(exp(c(-1.08, -1.08 + 1.48)), digits = 3),
           crazy = round(exp(c(-1.22, -1.22 + 1.27)), digits = 3),
           juice = round(exp(c(-1.21, -1.21 + 1.39)), digits = 3),
           chop = round(exp(c(-1.14, -1.14 + 1.41)), digits = 3))

# and look at the partially pooled mean
plot(hazard_model, pars = c("b_ret_mean",
                            "b_trt_ret_mean",
                            "b_pil_mean",
                            "b_trt_pil_mean"))

print(hazard_model, pars = c("b_ret_mean",
                             "b_trt_ret_mean",
                             "b_pil_mean",
                             "b_trt_pil_mean"))

# next: add hazard ratios and total effect (hazard ratio) for treatment effects
# in generated quantities
# AND tune priors!


#_______________________________________________________________________________________________
# 7. Save image ----
#_______________________________________________________________________________________________

save.image("RData - final/poisson_model.RData")
