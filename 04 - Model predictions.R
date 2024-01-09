# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04 - Model predictions
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 29 Dec 2023
# Date completed: 
# Date last modified: 09 Jan 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)
library(ggridges)        # ridgeline plot

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

load("poisson_model.RData")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(m5))

#_______________________________________________________________________________________________
# 4. Visualize time-varying baseline hazard function ----
#_______________________________________________________________________________________________
# 4a. Prepare predictions ----
#_______________________________________________________________________________________________

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
weeks.pred <- seq(1, 52, length.out = 100)           # 100 increments to predict on

basis.pred <- bs(weeks.pred,       
                 knots = knot.list[-c(1, n.knots)],          
                 degree = 3,                     
                 intercept = FALSE)

# we'll start by calculating predictions for every draw
spline.preds.df <- data.frame()

for (i in 1:nrow(model.draws)) {
  
  # subset spline parameters only
  spline.params <- model.draws[i , 1:7]
  
  # convert to vector for multiplication
  w0 <- t(as.matrix(as.numeric(spline.params[ ,2:7])))
  
  # multiply with 'sweep'
  w0.by.b <- sweep(basis.pred, 2, w0, `*`)
  
  # sum to create "normal" scale spline prediction
  w0.by.b.sum <- apply(w0.by.b, 1, sum)
  
  # add to the intercept (a0) and exponentiate to transform to hazard rate scale
  spline.pred <- exp(spline.params$a0 + w0.by.b.sum)
  
  # create df to hold predictions
  spline.pred.df <- data.frame(x = weeks.pred,
                               y = spline.pred,
                               draw = i)
  
  # bind into master df
  spline.preds.df <- bind_rows(spline.preds.df,
                               spline.pred.df)
  
}

# summarize the predictions with quantiles
spline.preds.df.med <- spline.preds.df %>%
  
  group_by(x) %>%
  
  summarise(q = list(quantile(y, probs = c(0.025,
                                           0.05,
                                           0.10,
                                           0.125,
                                           0.50,
                                           0.875,
                                           0.90,
                                           0.95,
                                           0.975)))) %>% 
  
  unnest_wider(q)

#_______________________________________________________________________________________________
# 4b. Plot median predictions with all draws ----
#_______________________________________________________________________________________________

ggplot() +
  
  # white background
  theme_bw() +
  
  # line for median spline prediction
  geom_line(data = spline.preds.df.med,
            aes(x = x,
                y = `50%`),     
            linewidth = 1.25) +
  
  # lines for all 4000 predictions
  geom_line(data = spline.preds.df,
            aes(x = x,
                y = y,
                group = draw),
            alpha = 0.01) +
  
  # axis labels
  ylab("Baseline hazard") +
  xlab("Week of year")

#_______________________________________________________________________________________________
# 4c. Plot median predictions with credible intervals ----
#_______________________________________________________________________________________________

# pivot longer for plotting with a legend
spline.preds.df.med.long <- spline.preds.df.med %>%
  
  pivot_longer(cols = 2:ncol(spline.preds.df.med)) %>%
  
  mutate(level = ifelse(name %in% c("2.5%", "97.5%"),
                        "95%",
                        ifelse(name %in% c("5%", "95%"),
                               "90%",
                               ifelse(name %in% c("12.5%", "87.5%"),
                                      "75%",
                                      "50%"))))

# add repeating pattern for lower and upper intervals
spline.pred.ci <- spline.preds.df.med.long %>% 
  
  filter(level != "50%") %>%
  
  mutate(bound = rep_len(c("lo", "lo", "lo", "up", "up", "up"), n())) %>%
  
  dplyr::select(-name) %>%
  
  pivot_wider(names_from = c("bound")) %>%
  
  mutate(level = factor(level, levels = c("95%", "90%", "75%")))

# plot
ggplot() +
  
  # white background
  theme_bw() +
  
  # lines for all 4000 predictions
  geom_ribbon(data = spline.pred.ci,
              aes(x = x,
                  ymin = lo,
                  ymax = up,
                  alpha = level)) +
  
  # define alpha levels
  scale_alpha_discrete(range = c(0.05, 0.25)) +
  
  # line for median spline prediction
  geom_line(data = spline.preds.df.med,
            aes(x = x,
                y = `50%`),     
            linewidth = 1.25) +
  
  # axis labels
  ylab("Baseline hazard") +
  xlab("Week of year")

#_______________________________________________________________________________________________
# 5. Visualize covariate effects and evaluate strength of evidence ----
#_______________________________________________________________________________________________
# 5a. All covariates (hazard ratios) - ridgeline plot ----
#_______________________________________________________________________________________________

# select parameters 
covariate.draws <- model.draws %>% dplyr::select(hr_sex:hr_bci)

# pivot
covariate.draws.long <- covariate.draws %>% 
  
  pivot_longer(cols = hr_sex:hr_bci)

# plot KDEs
ggplot(covariate.draws.long,
       aes(x = value,
           y = name)) +
  
  # white background
  theme_bw() +
  
  # KDE
  geom_density_ridges(fill = "lightgray") +
  
  # intercept at 1 (hazard ratio)
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  
  # coordinates
  coord_cartesian(xlim = c(0, 3.5)) +
  
  # axis labels
  xlab("Hazard ratio") +
  ylab("Parameter")

#_______________________________________________________________________________________________
# 5b. All covariates (hazard ratios) - point estimates [medians] and CIs ----
#_______________________________________________________________________________________________

# here we'll use 95% and 50% CIs
sex.ci <- data.frame(med = median(covariate.draws$hr_sex),
                     lo.1 = quantile(covariate.draws$hr_sex, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_sex, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_sex, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_sex, prob = 0.75),
                     var = "sex:M")

ret.ci <- data.frame(med = median(covariate.draws$hr_ret),
                     lo.1 = quantile(covariate.draws$hr_ret, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_ret, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_ret, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_ret, prob = 0.75),
                     var = "retention")

pil.ci <- data.frame(med = median(covariate.draws$hr_pil),
                     lo.1 = quantile(covariate.draws$hr_pil, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_pil, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_pil, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_pil, prob = 0.75),
                     var = "piling")

bsz.ci <- data.frame(med = median(covariate.draws$hr_bsz),
                     lo.1 = quantile(covariate.draws$hr_bsz, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_bsz, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_bsz, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_bsz, prob = 0.75),
                     var = "body size")

bci.ci <- data.frame(med = median(covariate.draws$hr_bci),
                     lo.1 = quantile(covariate.draws$hr_bci, prob = 0.025),
                     up.1 = quantile(covariate.draws$hr_bci, prob = 0.975),
                     lo.2 = quantile(covariate.draws$hr_bci, prob = 0.25),
                     up.2 = quantile(covariate.draws$hr_bci, prob = 0.75),
                     var = "body condition index")

# bind together
all.ci <- bind_rows(sex.ci, ret.ci, pil.ci, bsz.ci, bci.ci)

# plot
ggplot(data = all.ci,
       aes(x = med,
           y = var)) +
  
  # white background
  theme_bw() +
  
  # vertical line at 1
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  
  # credible intervals
  # 95%
  geom_errorbarh(aes(xmin = lo.1,
                     xmax = up.1,
                     y = var),
                 height = 0,
                 linewidth = 2,
                 color = "lightgray") +
  # 50%
  geom_errorbarh(aes(xmin = lo.2,
                     xmax = up.2,
                     y = var),
                 height = 0,
                 linewidth = 2,
                 color = "gray") +
  
  # points
  geom_point(aes(x = med,
                 y = var),
             size = 2.5) +
  
  # axis titles
  xlab("Hazard ratio") +
  ylab("") +
  
  # x-axis scale
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5))

#_______________________________________________________________________________________________
# 5c. All covariates (hazard ratios) - % of 95% interval on same side of 1 as median ----
#_______________________________________________________________________________________________

all.ci <- all.ci %>% 
  
  # total length of interval
  mutate(total.interval = up.1 - lo.1) %>%
  
  # which side of 1 is the median on?
  mutate(which.side = ifelse(med > 1,
                             1,
                             0)) %>%
  
  # determine what proportion is on the same side of 1
  mutate(v = ifelse(which.side == 1,
                    (up.1 - 1) / total.interval,
                    (1 - lo.1) / total.interval)) %>%
  
  # change every value > 1 
  mutate(v = ifelse(v > 1,
                    1,
                    v))

all.ci
