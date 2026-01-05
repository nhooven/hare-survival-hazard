# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 09 - Survivorship simulation
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 31 Dec 2025 
# Date completed: 05 Jan 2026
# Date last modified: 05 Jan 2026
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# We'll create a sample of simulated individuals
  # even sex ratio
  # average (with variability) body condition
  # different exposure to treatment

# we could start three treatments at once OR
# wait until some midway period to implement treatment

# at least to start with, we'll only use a subsample of posterior draws for efficiency
# we'll use model 2 because it's a nice balance

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)
library(tictoc)          # timing

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.2 <- read.csv("Model outputs/model_2.csv")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Body condition score distribution ----
#_______________________________________________________________________________________________

fates.byIndiv <- fates %>%
  
  group_by(deployment.1) %>%
  
  slice(1)

hist(fates.byIndiv$BCI.1)

possible.BCI <- (fates.byIndiv$BCI.1 - mean(fates$BCI.1)) / sd(fates$BCI.1)

hist(possible.BCI)

#_______________________________________________________________________________________________
# 4. Simulated individual datasets ----

# each forest type / treatment combination will include 100 animals
# we'll simulate for 2 years (104 weeks)

#_______________________________________________________________________________________________
# 4a. Full exposure ----

# treatments here start at week 1

#_______________________________________________________________________________________________

# SET 1a - SFL
set.1a <- data.frame(
  
  indiv = 1:300,
  trt = rep(c("control", "ret", "pil"), each = 100),
  sex = rep(c(0, 1), each = 50, times = 3),
  BCI.1 = sample(possible.BCI, size = 300)
  
) %>%
  
  # add sex_forest
  mutate(sex.forest = ifelse(sex == 0, 1, 2)) %>%
  
  # replicate rows
  slice(rep(1:n(), each = 104)) %>%
  
  # add week and scaled study.week
  mutate(
    
    week = rep(rep(1:52, times = 2), times = 300),
    
    study.week = (rep(1:104, times = 300) - mean(fates$study.week)) / sd(fates$study.week)
    
  ) %>%
  
  # add treatment variables
  mutate(
    
    post1 = case_when(
      
      week < 53 ~ 1,
      week > 52 ~ 0
      
    ),
    
    post2 = case_when(
      
      week > 52 ~ 1,
      week < 53 ~ 0
      
    ),
    
    ret = case_when(
      
      indiv %in% 1:100 ~ 0,
      indiv %in% 101:200 ~ 1,
      indiv %in% 201:300 ~ 0
      
    ),
    
    pil = case_when(
      
      indiv %in% 1:100 ~ 0,
      indiv %in% 101:200 ~ 0,
      indiv %in% 201:300 ~ 1
      
    )
    
  )

# SET 1b - XMC
set.1b <- data.frame(
  
  indiv = 1:300,
  trt = rep(c("control", "ret", "pil"), each = 100),
  sex = rep(c(0, 1), each = 50, times = 3),
  BCI.1 = sample(possible.BCI, size = 300)
  
) %>%
  
  # add sex_forest
  mutate(sex.forest = ifelse(sex == 0, 3, 4)) %>%
  
  # replicate rows
  slice(rep(1:n(), each = 104)) %>%
  
  # add week and scaled study.week
  mutate(
    
    week = rep(rep(1:52, times = 2), times = 300),
    
    study.week = (rep(1:104, times = 300) - mean(fates$study.week)) / sd(fates$study.week)
    
  ) %>%
  
  # add treatment variables
  mutate(
    
    post1 = case_when(
      
      week < 53 ~ 1,
      week > 52 ~ 0
      
    ),
    
    post2 = case_when(
      
      week > 52 ~ 1,
      week < 53 ~ 0
      
    ),
    
    ret = case_when(
      
      indiv %in% 1:100 ~ 0,
      indiv %in% 101:200 ~ 1,
      indiv %in% 201:300 ~ 0
      
    ),
    
    pil = case_when(
      
      indiv %in% 1:100 ~ 0,
      indiv %in% 101:200 ~ 0,
      indiv %in% 201:300 ~ 1
      
    )
    
  )

#_______________________________________________________________________________________________
# 4b. Partial exposure ----

# treatments here start at week 1
# but we start monitoring animals at week of year 40

# 01-05-2026 - We can try this later if people think it's necessary

#_______________________________________________________________________________________________

#_______________________________________________________________________________________________
# 5. Spline for prediction ----
#_______________________________________________________________________________________________

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

#_______________________________________________________________________________________________
# 6. Simulate lifetimes ----
#_______________________________________________________________________________________________
# 6a. Helper functions ----
#_______________________________________________________________________________________________

# kaplan meier curve
kap_mei <- function (
    
  x
  
) {
  
  # split by trt
  x.control <- x %>% filter(trt == "control")
  x.ret <- x %>% filter(trt == "ret")
  x.pil <- x %>% filter(trt == "pil")
  
  # generic function
  internal_kap_mei <- function (y) {
    
  # loop through all follow-up times
  all.S <- vector(length = max(y$t))
  s.i <- vector(length = max(y$t))
  all.ci.low <- vector(length = max(y$t))
  all.ci.upp <- vector(length = max(y$t))
  
  all.S[1] = 1.0 # initialize at 100%
  s.i[1] = 1.0 # initialize at 100%
  all.ci.low[1] <- 1.0
  all.ci.upp[1] <- 1.0
  
  for (i in 2:max(y$t)) {
    
    # subset
    indivs.i <- y[y$t == i, ]
    
    # calculate survival at interval
    s.i[i] = 1 - (sum(indivs.i$event) / nrow(indivs.i))
    
    # cumulative survival at interval
    all.S[i] = prod(s.i[1:i])
    
    # calculate 90 % confidence interval
    focal.var = ((all.S[i]^2) * (1 - all.S[i])) / nrow(indivs.i)
    
    all.ci.low[i] = all.S[i] - (1.64 * sqrt(focal.var))
    all.ci.upp[i] = all.S[i] + (1.64 * sqrt(focal.var))
    
  }
  
  # df for returning
  km.df <- data.frame(trt = y$trt[1],
                      t = 1:max(y$t),
                      S = all.S,
                      ci.low = all.ci.low,
                      ci.upp = all.ci.upp)
  
  return(km.df)
    
  }
  
  km.df.all <- rbind(internal_kap_mei(x.control),
                     internal_kap_mei(x.ret),
                     internal_kap_mei(x.pil))
  
  return(km.df.all)
  
}

# calculate Poisson intensity function
calc_intensity_2 <- function (
    
  y,    # covariates
  v     # coefficients
  
) {
  
  # BLH spline prediction
  # weights w
  w <- t(as.matrix(as.numeric(c(v[paste0("w", ".", y$sex.forest[1], "..1.")],
                                v[paste0("w", ".", y$sex.forest[1], "..2.")],
                                v[paste0("w", ".", y$sex.forest[1], "..3.")],
                                v[paste0("w", ".", y$sex.forest[1], "..4.")],
                                v[paste0("w", ".", y$sex.forest[1], "..5.")],
                                v[paste0("w", ".", y$sex.forest[1], "..6.")],
                                v[paste0("w", ".", y$sex.forest[1], "..7.")],
                                v[paste0("w", ".", y$sex.forest[1], "..8.")],
                                v[paste0("w", ".", y$sex.forest[1], "..9.")]))))
  
  # multiply and sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
  
  # add to intercept and exponentiate
  blh = as.numeric(exp(v[paste0("a0.", y$sex.forest, ".")] + w.by.b.sum[y$week]))
  
  # coefficient prediction (must be on the log scale!)
  coef.pred = exp(log(v$hr_bci) * y$BCI.1 +
                    log(v$hr_bci_study_week) * y$study.week * y$BCI.1 +
                    log(v$hr_pil_total1) * y$post1 * y$pil +
                    log(v$hr_pil_total2) * y$post2 * y$pil +
                    log(v$hr_ret_total1) * y$post1 * y$ret +
                    log(v$hr_ret_total2) * y$post2 * y$ret)
  
  # total intensity
  intens = blh * coef.pred
  
  # Poisson draw
  pois.draw = rpois(n = 1, intens)
  
  # pack into a df
  y.intens.pois <- data.frame(intens = intens,
                              pois.draw = rpois(length(intens), intens))
  
  return(y.intens.pois)
  
}

# function to simulate lifetimes
sim_lifetime_byDeploy <- function (z, v, i) {
  
  # calculate intensity and bind
  z.1 <- cbind(z, calc_intensity_2(y = z,
                                   v))
  
  # extract lifetime
  focal.lifetime <- data.frame(
    
    iter = i,
    indiv = z.1$indiv[1],
    trt = z.1$trt[1],
    lifetime = ifelse(
      
      1 %in% z.1$pois.draw,
      which(z.1$pois.draw > 0)[1],
      104                                 
      
    )
    
  )
  
  # return
  return(focal.lifetime)
  
}

#_______________________________________________________________________________________________
# 6b. Loop by iteration ----

# split by indiv
set.1a.split <- split(set.1a, set.1a$indiv)
set.1b.split <- split(set.1b, set.1b$indiv)

draws = model.fit.2 %>% slice(sample(1:n(), size = 1000))

#_______________________________________________________________________________________________

# function
sim_survivorship <- function (
    
  draws,
  set.split
  
) {

# loop through draws
all.km.curves <- data.frame()

for (i in 1:nrow(draws)) {
  
  iter.draw <- draws[i, ]
  
  # apply function to each deployment
  all.iter.lifetimes <- do.call(rbind, lapply(set.split, 
                                              sim_lifetime_byDeploy,
                                              v = iter.draw,
                                              i = i))
  
  
  # calculate K-M curve
  # create dataset
  
  # blank df
  km.data <- data.frame()
  
  # loop
  for (j in 1:nrow(all.iter.lifetimes)) {
    
    # vector of zero observations + 1
    zero.vector <- vector(length = all.iter.lifetimes$lifetime[all.iter.lifetimes$indiv == j] - 1)
    zero.vector[1:length(zero.vector)] <- 0
    zero.vector[length(zero.vector) + 1] <- 1
    
    # pack into df
    focal.df <- data.frame(
      
      j = j,
      trt = all.iter.lifetimes$trt[j],
      t = 1:length(zero.vector),
      event = zero.vector
      
    )
    
    km.data <- rbind(km.data, focal.df)
    
  }
  
  # apply KM function
  km.curve.data <- kap_mei(x = km.data)
  
  # add in iteration
  km.curve.data$iter = i
  
  # pack into df
  all.km.curves <- rbind(all.km.curves, km.curve.data)
  
}
  
  # return
  return(all.km.curves)

}

# apply function
tic()
set.1a.km <- sim_survivorship(draws, set.1a.split)
toc()

tic()
set.1b.km <- sim_survivorship(draws, set.1b.split)
toc()

# ~ 22 seconds for 10 iterations
# 1000 iterations would take ~ 40 minutes

#_______________________________________________________________________________________________
# 7. Summarize survivorship ----
#_______________________________________________________________________________________________
# 7a. Helper functions ----
#_______________________________________________________________________________________________

# extract survival
extract_surv <- function (
    
  set.km,
  weeks = c(26, 52, 78),
  ci = 0.95) {
  
  set.km.summary <- set.km %>%
    
    # filter by weeks
    filter(t %in% weeks) %>%
    
    group_by(t, trt) %>%
    
    # summarize
    summarize(
      
      median = median(S),
      
      low = quantile(S, probs = (1 - ci) / 2),
      
      upp = quantile(S, probs = 1 - (1 - ci) / 2),
      
    )
  
  return(set.km.summary)
  
}

# prediction envelopes
km_pred <- function (x) {
  
  x.1 <- x %>%
    
    group_by(trt, t) %>%
    
    mutate(med = median(S),
           pe.low = quantile(S, prob = 0.025),
           pe.upp = quantile(S, prob = 0.975)) %>%
    
    slice(1) %>%
    
    ungroup() %>%
    
    dplyr::select(trt, t, med, pe.low, pe.upp)
  
  return(x.1)
  
}

#_______________________________________________________________________________________________
# 7b. Summarize and plot CIs ----
#_______________________________________________________________________________________________

set.1.summary <- rbind(
  
  extract_surv(set.1a.km) %>%
    
    mutate(forest = "SFL"),
  
  extract_surv(set.1b.km) %>%
    
    mutate(forest = "XMC")
  
)

# labels
set.1.summary <- set.1.summary %>%
  
  mutate(t = factor(t, levels = c(26, 52, 78),
                    labels = c("26 weeks", 
                                  "52 weeks",
                                  "78 weeks")),
         forest = factor(forest, levels = c("SFL",
                                            "XMC"),
                                 labels = c("spruce-fir-lodgepole",
                                            "xeric mixed conifer")))

# plot
ggplot(data = set.1.summary) +
  
  theme_bw() +
  
  facet_grid(forest ~ t) +
  
  geom_errorbar(aes(x = trt,
                    ymin = low,
                    ymax = upp,
                    color = trt),
                width = 0,
                linewidth = 1.5) +
  
  geom_point(aes(x = trt,
                 y = median,
                 fill = trt),
             shape = 21,
             size = 3) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0)) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  ylab("Cumulative survival") +
  
  scale_x_discrete(labels = c("control", "piling", "retention")) +
  
  scale_color_manual(values = c("gray75", "orange", "purple")) +
  scale_fill_manual(values = c("gray75", "orange", "purple"))

#_______________________________________________________________________________________________
# 7c. Plot KM curves ----
#_______________________________________________________________________________________________

# plot set 1 together
set.1a.km.pred <- km_pred(set.1a.km)
set.1b.km.pred <- km_pred(set.1b.km)

# add forest type
set.1a.km$forest <- "SFL"
set.1b.km$forest <- "XMC"

set.1a.km.pred$forest <- "SFL"
set.1b.km.pred$forest <- "XMC"

# bind together
set.1.km <- rbind(set.1a.km, set.1b.km)
set.1.km.pred <- rbind(set.1a.km.pred, set.1b.km.pred)

# plot
ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ forest,
             nrow = 2) +
  
  # prediction intervals
  geom_ribbon(data = set.1.km.pred,
              aes(x = t,
                  y = med,
                  ymin = pe.low,
                  ymax = pe.upp,
                  group = trt,
                  fill = trt),
              alpha = 0.25) +
  
  geom_line(data = set.1.km.pred,
            aes(x = t,
                y = med,
                linetype = trt,
                color = trt),
            linewidth = 0.85) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.background = element_rect(color = "black"),
        legend.position = c(0.85, 0.30),
        legend.title = element_blank(),
        # text inside panels
        strip.text = element_text(hjust = 0.98,
                                  margin = margin(b = -15, 
                                                  l = 10,
                                                  t = 0),
                                  size = 14),
        strip.background = element_blank(),
        strip.clip = "off") +
  
  xlab("Week") +
  ylab("Cumulative survival") +
  
  coord_cartesian(xlim = c(5.5, 97)) +
  
  scale_x_continuous(breaks = c(13, 26, 39, 52, 65, 78, 92)) +
  
  scale_color_manual(values = c("gray75", "orange", "purple")) +
  scale_fill_manual(values = c("gray75", "orange", "purple"))

# well, doesn't look like treatments do much here
# next, we could start all our animals at year 2
