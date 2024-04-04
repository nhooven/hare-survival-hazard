# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 04c - Model predictions (survival functions)
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 11 Mar 2024
# Date completed: 11 Mar 2024
# Date last modified: 11 Mar 2024
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(rstan)           # modeling with Stan
library(splines)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

load("RData - final/poisson_model.RData")

#_______________________________________________________________________________________________
# 3. Extract parameter estimates ----
#_______________________________________________________________________________________________

# extract draws
model.draws <- as.data.frame(rstan::extract(m7))

#_______________________________________________________________________________________________
# 4. Calculate cumulative survival curves ----
#_______________________________________________________________________________________________
# 4a. Prepare predictions ----
#_______________________________________________________________________________________________

# wekk-by-week predictions
# use one year for prediction
weeks.pred <- seq(1, 52, length.out = 100)           # 100 increments to predict on

# we'll begin this ~ October 9 (week 41)
weeks.pred.1 <- weeks.pred + 40 
weeks.pred.1 <- ifelse(weeks.pred.1 > 52,
                       weeks.pred.1 - 52,
                       weeks.pred.1)

weeks.pred.1

# basis functions
basis.pred <- bs(weeks.pred.1,       
                 knots = knot.list[-c(1, n.knots)],          
                 degree = 3,                     
                 intercept = FALSE)

# we'll calculate predictions for every draw
preds.df <- data.frame()

for (i in 1:nrow(model.draws)) {
  
  # subset model draws
  focal.draws <- model.draws[i, ]
  
  # subset spline parameters only
  spline.params <- focal.draws[ , c(1, 7:12)]
  
  # convert to vector for multiplication
  w0 <- t(as.matrix(as.numeric(spline.params[ , 2:7])))
  
  # multiply with 'sweep'
  w0.by.b <- sweep(basis.pred, 2, w0, `*`)
  
  # sum to create "normal" scale spline prediction
  w0.by.b.sum <- apply(w0.by.b, 1, sum)
  
  # add to the intercept (a0) and exponentiate to transform to hazard rate scale
  spline.pred <- exp(spline.params$a0 + w0.by.b.sum)
  
  # and now the covariate effects
  cov.effects.ret <- exp(focal.draws$b_sex * 0 +
                         focal.draws$b_mas * 0 +
                         focal.draws$b_hfl * 0 +
                         focal.draws$b_bci * 0 +
                         focal.draws$b_ret * 1 +
                         focal.draws$b_pil * 0 +
                         focal.draws$b_sex * 0)
  
  cov.effects.pil <- exp(focal.draws$b_sex * 0 +
                         focal.draws$b_mas * 0 +
                         focal.draws$b_hfl * 0 +
                         focal.draws$b_bci * 0 +
                         focal.draws$b_ret * 0 +
                         focal.draws$b_pil * 1 +
                         focal.draws$b_sex * 0)
  
  # calculate hazard, cumulative hazard, and cumulative survival
  # weekly hazard
  hazard.control <- spline.pred
  hazard.ret <- spline.pred * cov.effects.ret
  hazard.pil <- spline.pred * cov.effects.pil
  
  # cumulative hazard
  hazard.control.cumul <- cumsum(spline.pred)
  hazard.ret.cumul <- cumsum(hazard.ret)
  hazard.pil.cumul <- cumsum(hazard.pil)
  
  # cumulative survival
  hazard.control.surv <- exp(-hazard.control.cumul)
  hazard.ret.surv <- exp(-hazard.ret.cumul)
  hazard.pil.surv <- exp(-hazard.pil.cumul)
  
  # create df to hold predictions
  pred.df <- data.frame(x = weeks.pred.1,
                        x.1 = weeks.pred,     # ordered correctly
                        control = hazard.control.surv,
                        ret = hazard.ret.surv,
                        pil = hazard.pil.surv,
                        draw = i)
  
  # bind into master df
  preds.df <- bind_rows(preds.df, pred.df)
  
}

# pivot longer
preds.df.long <- preds.df %>% pivot_longer(cols = 3:5)

# summarize the predictions with quantiles
preds.df.med <- preds.df.long %>%
  
  group_by(x.1, name) %>%
  
  summarise(median = median(value),
            q = list(quantile(value, probs = c(0.025, 0.975)))) %>% 
  
  unnest_wider(q)

#_______________________________________________________________________________________________
# 5. Plot median predictions (with credible intervals) ----
#_______________________________________________________________________________________________

# label factor levels
preds.df.med$name <- factor(preds.df.med$name,
                            labels = c("control",
                                       "piling",
                                       "retention"))

ggplot(data = preds.df.med) +
  
  theme_bw() +
  
  # credible intervals
  geom_ribbon(aes(x = x.1,
                  ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = name),
              alpha = 0.15) +
  
  # median predictions
  geom_line(aes(x = x.1,
                y = median,
                color = name,
                linetype = name),
            linewidth = 1.25) +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.7),
        legend.title = element_blank()) +
  
  # axis labels
  ylab("Survival") +
  xlab("Week of the year") +
  
  # coordinates
  coord_cartesian(xlim = c(3.5, 49.7)) +
  
  # scale x axis correctly
  scale_x_continuous(breaks = c(2, 12, 22, 32, 42, 52),
                     labels = c(43, 1, 11, 21, 31, 41)) +
  
  # colors
  scale_color_manual(values = c("black", "purple", "orange")) +
  scale_fill_manual(values = c("black", "purple", "orange")) +
  
  # linetypes 
  scale_linetype_manual(values = c("solid", "dotted", "dashed"))
