# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 10 - Cumulative survival
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 01 May 2026
# Date completed: 01 May 2026
# Date last modified: 01 May 2026
# R version: 4.4.3

#_______________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________

# instead of simulations (which make sense for PPCs), let's just calculate the
# implied cumulative survival for each treatment,
# beginning Oct 1 and looking at 3-month, 6-month, and 12-month survival probability

# we'll do year 1 and year 2 post separately

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)            # splines
library(bayestestR)      # HDIs

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

# model samples
model.fit.2 <- readRDS("models/model_2.rds")

# dataset
fates <- read.csv("data_cleaned/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("data_cleaned/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("data_cleaned/day_lookup.csv")

#_______________________________________________________________________________
# 3. Initialize data ----

# these can all go in the same data.frame
# pattern is:
#   F: post1 C, post2 C, post1 R, post2 R, post1 P, post2 P
#   M: post1 C, post2 C, post1 R, post2 R, post1 P, post2 P

#_______________________________________________________________________________

set.to.pred <- data.frame(
  
  sex = c(0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1),
  post1 = c(1, 0, 1, 0, 1, 0,
            1, 0, 1, 0, 1, 0),
  post2 = c(0, 1, 0, 1, 0, 1,
            0, 1, 0, 1, 0, 1),
  ret = c(0, 0, 1, 1, 0, 0,
          0, 0, 1, 1, 0, 0),
  pil = c(0, 0, 0, 0, 1, 1,
          0, 0, 0, 0, 1, 1),
  
  # keep at their means
  BCI.1 = 0,
  DM = 0,
  OPEN = 0,
  
  # index
  index = 1:12
  
)

#_______________________________________________________________________________
# 4. Spline for prediction ----
#_______________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# weeks to predict on
weeks.pred <- 1:52

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________
# 5. Calculate cumulative survival - function ----

# x is each covariate configuration from set.to.pred
# y is the current iteration

#_______________________________________________________________________________

calc_surv <- function (x, y) {
  
  # expand x to include all weeks in a year
  suppressWarnings(
    
    x.1 <- cbind(x, week = 1:52)
    
  )
  
  # 1. calculate hazard for every week
  # extract spline weights - ws[s, b] -  just the means 
  w <- t(
    
    as.matrix(
      
      as.numeric(
        
        c(
          
          y[paste0("ws", "[", x$sex[1] + 1, ", 1]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 2]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 3]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 4]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 5]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 6]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 7]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 8]")],
          y[paste0("ws", "[", x$sex[1] + 1, ", 9]")]
          
        )
        
      )
      
    )
    
  )
  
  # multiply and sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
  
  # add to intercept and exponentiate
  blh = exp(
    
    as.numeric(
      
      y[paste0("a0s[", x$sex[1] + 1, "]")]
      
    ) + w.by.b.sum
    
  )
  
  # total hazard ratio prediction
  hr = exp(
    
    # all should be zero
    y["b_bci"] * x$BCI +
    y["b_dm"] * x$DM +
    y["b_open"] * x$OPEN +
      
    y["b_post1"] * x$post1 +
    y["b_post2"] * x$post2 +
    y["b_ret_post1"] * x$ret * x$post1 +
    y["b_ret_post2"] * x$ret * x$post2 +
    y["b_pil_post1"] * x$pil * x$post1 +
    y["b_pil_post2"] * x$pil * x$post2
    
    )
  
  # full hazard
  full.haz = blh * hr
  
  # sum hazard for each period, add to data.frame
  x <- x %>%
    
    mutate(
      
      month.3 = exp(-sum(full.haz[1:(3 * 4)])),
      month.6 = exp(-sum(full.haz[1:(6 * 4)])),
      month.12 = exp(-sum(full.haz[1:(12 * 4)]))
      
    )
  
  # returns a 1-row df
  return(x)
  
}

#_______________________________________________________________________________
# 6. Calculate cumulative survival - loop ----

set.to.pred.split <- split(set.to.pred, set.to.pred$index)

draws <- do.call(rbind, model.fit.2)

#_______________________________________________________________________________

all.S <- data.frame()

for (i in 1:nrow(draws)) {
  
  iter.draw <- draws[i, ]
  
  # apply function to all splits
  set.to.pred.focal <- lapply(set.to.pred.split, 
                              calc_surv,
                              y = iter.draw) |> 
    
    # bind
    do.call(rbind, args = _)
  
  # bind to master
  all.S <- rbind(all.S, set.to.pred.focal)
  
}

# summarize with posterior medians, CIs, and prepare for plotting
all.S.summary <- all.S %>%
  
  group_by(index) %>%
  
  summarize(
    
    # medians
    med.3 = median(month.3),
    med.6 = median(month.6),
    med.12 = median(month.12),
    
    # l90
    l90.3 = as.numeric(hdi(month.3, ci = 0.90)[2]),
    l90.6 = as.numeric(hdi(month.6, ci = 0.90)[2]),
    l90.12 = as.numeric(hdi(month.12, ci = 0.90)[2]),
    
    # u90
    u90.3 = as.numeric(hdi(month.3, ci = 0.90)[3]),
    u90.6 = as.numeric(hdi(month.6, ci = 0.90)[3]),
    u90.12 = as.numeric(hdi(month.12, ci = 0.90)[3])
    
  ) %>%
  
  # join in original
  left_join(
    
    set.to.pred,
    by = "index"
    
  ) %>%
  
  # factors
  mutate(
    
    sex = ifelse(sex == 0, "F", "M"),
    yr = case_when(post1 == 1 ~ "YR 1",
                   post2 == 1 ~ "YR 2"),
    trt = case_when(ret == 0 & pil == 0 ~ "control",
                    ret == 1 & pil == 0 ~ "retention",
                    ret == 0 & pil == 1 ~ "piling")
    
  ) %>%
  
  # drop unused columns
  dplyr::select(-c(post1, post2, ret, pil, BCI.1, DM, OPEN, index)) %>%
  
  # pivot properly
  pivot_longer(
    
    cols = med.3:u90.12,
    names_to = c("var", "month"),
    names_sep = 4
    
  ) %>% 
  
  pivot_wider(
    
    names_from = var
    
  ) %>%
  
  # month as factor
  mutate(month = factor(month,
                        levels = c("3", "6", "12")))

#_______________________________________________________________________________
# 7. Plot ----
#_______________________________________________________________________________

ggplot(data = all.S.summary) +
  
  theme_bw() +
  
  facet_grid(yr ~ trt) +
  
  geom_errorbar(aes(x = month,
                    y = med.,
                    ymin = l90.,
                    ymax = u90.,
                    color = sex),
                width = 0,
                position = position_dodge(width = 0.5),
                linewidth = 1.1,
                alpha = 0.60) +
  
  geom_point(aes(x = month,
                 y = med.,
                 fill = sex),
             position = position_dodge(width = 0.5),
             shape = 21) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        legend.position = c(0.76, 0.65),
        legend.title = element_blank(),
        legend.background = element_rect(color = NA,
                                         fill = NA)) +
  
  # labels
  xlab("Months") +
  ylab("Survival probability") +
  
  # colors
  scale_color_manual(values = c("#FF3300", "gray25")) +
  scale_fill_manual(values = c("#FF3300", "gray25"))

# 512 x 353

#_______________________________________________________________________________
# 8. Summary table ----
#_______________________________________________________________________________

write.table(all.S.summary, "clipboard", sep = "\t")
