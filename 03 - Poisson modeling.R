# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03b - Poisson modeling
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 03 Dec 2023
# Date completed: 29 Dec 2023
# Date last modified: 06 Oct 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)         # construct basic functions
library(mgcv)            # cyclic splines

#_______________________________________________________________________________________________
# 2. Read in and format data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_cleaned_09_28_2024.csv")

# keep only columns we need for modeling
fates.1 <- fates %>% dplyr::select(cluster,
                                   Site,
                                   Ear.tag,
                                   week,
                                   year,
                                   mort,
                                   cens,
                                   Collar.type.1,
                                   Sex.1,
                                   Mass.1,
                                   HFL.1,
                                   PC1,
                                   BCI.1,
                                   post.trt,
                                   trt.ret,
                                   trt.pil)

#_______________________________________________________________________________________________
# 3. Construct basis functions for spline on hazard ----
#_______________________________________________________________________________________________

# procedure adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

# define knots (quantile)
n.knots <- 5 + 1

knot.list <- quantile(fates.1$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# creating the basis matrix (cyclic spline)
basis <- cSplineDes(fates.1$week,
                    knots = knot.list,
                    ord = 4)                # cubic

# plot basis functions
basis.plot <- tibble(week = fates.1$week,
                     b1 = basis[ , 1],
                     b2 = basis[ , 2],
                     b3 = basis[ , 3],
                     b4 = basis[ , 4],
                     b5 = basis[ , 5]) %>%
              pivot_longer(cols = b1:b5)

ggplot(data = fates.1,
       aes(x = week,
           y = mort)) +
  theme_bw() +
  geom_jitter(alpha = 0.1,
              height = 0.01) +
  geom_line(data = basis.plot,
            aes(x = week,
                y = value,
                group = name),
            linewidth = 1.5,
            alpha = 0.25)

#_______________________________________________________________________________________________
# 4. Build data list ----
#_______________________________________________________________________________________________

fates.stan.1 <- list(N = nrow(fates.1),
                     y_mort = fates.1$mort,
                     y_cens = fates.1$cens,
                     t = fates.1$week,
                     collar = fates.1$Collar.type.1,
                     sex = fates.1$Sex.1,
                     mas = fates.1$Mass.1,
                     hfl = fates.1$HFL.1,
                     bsz = fates.1$PC1,
                     bci = fates.1$BCI,
                     trt = fates.1$post.trt,
                     ret = fates.1$trt.ret,
                     pil = fates.1$trt.pil,
                     clust = fates.1$cluster,
                     nclust = 4,
                     basis = basis,
                     nbasis = ncol(basis))

#_______________________________________________________________________________________________
# 5. Run first model (random intercept for cluster) ----
#_______________________________________________________________________________________________

m1 <- rstan::stan(
  
  file = "m1.stan",
  data = fates.stan.1,
  chains = 1,
  warmup = 1000,
  iter = 2000
  
)

print(m1)

# estimates
plot(m1, pars = c("b_sex", "b_ret", "b_pil", "b_trt_r", "b_trt_p", "b_mas", "b_hfl", "b_bci"))
plot(m1, pars = c("a_c1", "a_c2", "a_c3", "a_c4"))

# trace
traceplot(m1, pars = c("b_sex", "b_ret", "b_pil", "b_trt_r", "b_trt_p", "b_mas", "b_hfl", "b_bci"))
traceplot(m1, pars = c("a_c1", "a_c2", "a_c3", "a_c4"))

#_______________________________________________________________________________________________
# 6. Run second model (censoring likelihood) ----
#_______________________________________________________________________________________________

m2 <- rstan::stan(
  
  file = "m2.stan",
  data = fates.stan.1,
  chains = 4,
  warmup = 1000,
  iter = 2000
  
)

print(m2)

# estimates
plot(m2, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_mas", "hr_hfl", "hr_bci"))
plot(m2, pars = c("hr_col"))

# trace
traceplot(m2, pars = c("hr_sex", "hr_ret", "hr_pil", "hr_mas", "hr_hfl", "hr_bci"))
traceplot(m2, pars = c("cens0", "hr_col"))

#_______________________________________________________________________________________________
# 7. Run third model (censoring likelihood, constraining spline weights even more) ----
#_______________________________________________________________________________________________

# here we need to balance (1) the seasonal variance in hazard accounted for by
# the spline and (2) the variance attributable to covariates

# initial N(0, 1) priors for spline weights were very uncertain and led to 
# a large spike in spring mortality
# balancing the magnitude and spread of these coefficients will be important
# for inference - we don't want to overfit! 
# more skeptical priors should draw our spline function closer to zero

m3 <- rstan::stan(
  
  file = "m3.stan",
  data = fates.stan.1,
  chains = 1,
  warmup = 1000,
  iter = 2000
  
)

print(m3)

# estimates
plot(m3, pars = c("w0"))
plot(m3, pars = c("b_sex", "b_ret", "b_pil", "b_mas", "b_hfl", "b_bci"))
plot(m3, pars = c("b_col"))


# for these interactive effects, calculate the "total coefficient" as:
# (b_ret * ret) + (b_trt_r * trt * ret), etc. 

# show these effects as "before-after" plots

#_______________________________________________________________________________________________
# 7. Save image ----
#_______________________________________________________________________________________________

save.image("RData - final/poisson_model.RData")
