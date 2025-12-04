# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06b - Posterior predictive checks - survivorship curves
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 04 Dec 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# Simulated lifetimes
# Here we'll simulate lifetimes from a subset of posterior draws, create Kaplan-Meier
# survivorship curves for each "population", then compare to the empirical K-M curve
# I'll just need to snap every individual to the same starting place

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
model.fit.1 <- read.csv("Model outputs/model_1.csv")
model.fit.2 <- read.csv("Model outputs/model_2.csv")
model.fit.3 <- read.csv("Model outputs/model_3.csv")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# day lookup table
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________________________
# 3a. Attribute study-week to fates df ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  # join in year
  mutate(year = fates.year$year) %>%
  
  left_join(
    
    day.lookup %>%
      
      dplyr::select(
        
        year,
        study.year.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1),
         study.week.s = (study.week - mean(study.week)) / sd(study.week))

#_______________________________________________________________________________________________
# 3b. Spline for prediction ----
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
# 4. Split dataset by deployment ----
#_______________________________________________________________________________________________

# split by deployment
fates.1.deploy.split <- split(fates.1, fates.1$deployment.1)

#_______________________________________________________________________________________________
# 5. Prepare dataset for Kaplan-Meier curves ----
#_______________________________________________________________________________________________
# 5a. Prepare empirical data ----

# Here we'll keep everything the same, but include a "t" variable denoting the week from the beginning of 
# monitoring per deployment

#_______________________________________________________________________________________________

# function to incorporate t since monitoring started
add_t <- function (x) {
  
  x$t = 1:nrow(x)
  
  return(x)
  
}

# apply
fates.2 <- do.call(rbind, lapply(fates.1.deploy.split, add_t))

#_______________________________________________________________________________________________
# 5b. Write function ----

# Ideally we can also use this within the simulation procedure
# we'll need to split data based on t since deployment

#_______________________________________________________________________________________________

kap_mei <- function (
    
  x,
  response = "y.mort.scen1"
  
) {
  
  # loop through all follow-up times
  all.S <- vector(length = max(x$t))
  s.i <- vector(length = max(x$t))
  all.ci.low <-vector(length = max(x$t))
  all.ci.upp <-vector(length = max(x$t))
  
  all.S[1] = 1.0 # initialize at 100%
  s.i[1] = 1.0 # initialize at 100%
  all.ci.low[1] <- 1.0
  all.ci.upp[1] <- 1.0
  
  for (i in 2:max(x$t)) {
    
    # subset
    indivs.i <- x[x$t == i, ]
    
    # calculate survival at interval
    s.i[i] = 1 - (sum(indivs.i[[response]]) / nrow(indivs.i))
    
    # cumulative survival at interval
    all.S[i] = prod(s.i[1:i])
    
    # calculate 90 % confidence interval
    focal.var = ((all.S[i]^2) * (1 - all.S[i])) / nrow(indivs.i)
    
    all.ci.low[i] = all.S[i] - (1.64 * sqrt(focal.var))
    all.ci.upp[i] = all.S[i] + (1.64 * sqrt(focal.var))
    
  }
  
  # df for returning
  km.df <- data.frame(t = 1:max(x$t),
                      S = all.S,
                      ci.low = all.ci.low,
                      ci.upp = all.ci.upp)
  
  return(km.df)
  
}

# apply function
# model 1
km.df.1 <- kap_mei(fates.2,
                   response = "y.mort.scen1")

# model 2
km.df.2 <- kap_mei(fates.2,
                   response = "y.mort.scen2")

# model 3
km.df.3 <- kap_mei(fates.2,
                   response = "y.mort.scen3")

# plot each
ggplot() +
  
  theme_bw() +
  
  # scen1
  geom_ribbon(data = km.df.1,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "gray") +
  
  geom_line(data = km.df.1,
            aes(x = t,
                y = S)) +
  
  # scen2
  geom_ribbon(data = km.df.2,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "red") +
  
  geom_line(data = km.df.2,
            aes(x = t,
                y = S),
            color = "darkred") +
  
  # scen3
  geom_ribbon(data = km.df.3,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25,
              fill = "gold") +
  
  geom_line(data = km.df.3,
            aes(x = t,
                y = S),
            color = "orange") +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1))

#_______________________________________________________________________________________________
# 6. Simulate lifetimes ----

# Here we'll loop through each iteration (take a subsample first),
# loop through each deployment (n = 383),
# applying a calc intensity function on an 86-week dataset,
# take a Poisson draw and truncate the dataset after the first event
# then apply the K-M curve function

# the output will be a collection of K-M curves

#_______________________________________________________________________________________________
# 6a. Helper functions ----
#_______________________________________________________________________________________________

# calculate Poisson intensity function
calc_intensity_2 <- function (
    
  y,    # covariates
  v     # coefficients
  
  ) {
  
  # BLH spline prediction
  # weights w
  w <- t(as.matrix(as.numeric(c(v[paste0("w", ".", y$sex_forest[1], "..1.")],
                                v[paste0("w", ".", y$sex_forest[1], "..2.")],
                                v[paste0("w", ".", y$sex_forest[1], "..3.")],
                                v[paste0("w", ".", y$sex_forest[1], "..4.")],
                                v[paste0("w", ".", y$sex_forest[1], "..5.")],
                                v[paste0("w", ".", y$sex_forest[1], "..6.")],
                                v[paste0("w", ".", y$sex_forest[1], "..7.")],
                                v[paste0("w", ".", y$sex_forest[1], "..8.")],
                                v[paste0("w", ".", y$sex_forest[1], "..9.")]))))
  
  # multiply and sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
  
  # add to intercept and exponentiate
  blh = as.numeric(exp(v[paste0("a0.", y$sex_forest, ".")] + w.by.b.sum[y$week]))
  
  # coefficient prediction (must be on the log scale!)
  coef.pred = exp(log(v$hr_bci) * y$BCI.s +
                  log(v$hr_bci_study_week) * y$study.week.s * y$BCI.s +
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
sim_lifetime_byDeploy <- function (z) {
  
  # create new data.frame to fill
  suppressMessages(
    
    deploy.df <- data.frame(
      
      deployment.1 = z$deployment.1[1],
      sex_forest = z$sex_forest[1],
      t = 1:200,
      study.week = seq(z$study.week[1], z$study.week[1] + 19),
      post1 = z$post1[1],
      post2 = z$post2[1],
      ret = z$ret[1],
      pil = z$pil[1],
      BCI.s = z$BCI.s[1]
      
    ) %>%
      
      mutate(
        
        # subtract 160 (max week) from anything above it in study.week to recycle
        study.week = ifelse(study.week > 160,
                            (study.week - 160) + 4,    # so the correct week of the year will be attributed
                            study.week)
        
      ) %>%
      
      mutate(
        
        # standardize study.week
        study.week.s = (study.week - mean(fates.1$study.week)) / sd(fates.1$study.week)
        
      ) %>%
      
      # attribute correct week of the year
      left_join(
        
        day.lookup %>%
          
          dplyr::select(study.year.week,
                        study.week) %>%
          
          rename(week = study.year.week) %>%
          
          group_by(study.week) %>%
          
          slice(1)
        
      )
    
  )
  
  # calculate intensity and bind
  deploy.df <- cbind(deploy.df, calc_intensity_2(y = deploy.df,
                                                 v = iter.draw))
  
  # extract lifetime
  focal.lifetime <- data.frame(
    
    iter = i,
    deployment = z$deployment.1[1],
    lifetime = ifelse(
      
      1 %in% deploy.df$pois.draw,
      which(deploy.df$pois.draw > 0)[1],
      200                                 # we'll probably truncate this dataset
      # since the max weeks we monitored a hare was 86
    )
    
  )
  
  # return
  return(focal.lifetime)
  
}

#_______________________________________________________________________________________________
# 6b. Loop by iteration ----

# split deployments
deploy.split <- split(fates.2, fates.2$deployment.1)

draws = model.fit.1

#_______________________________________________________________________________________________

# loop through draws
all.km.curves <- data.frame()

for (i in 1:nrow(draws)) {
  
  iter.draw = draws[i, ]
  
  # apply function to each deployment
  all.iter.lifetimes <- do.call(rbind, lapply(deploy.split, 
                                              sim_lifetime_byDeploy))
  
  
  # calculate K-M curve
  # create dataset
  
  # blank df
  km.data <- data.frame()
  
  # loop
  for (j in 1:nrow(all.iter.lifetimes)) {
    
    # vector of zero observations + 1
    zero.vector <- vector(length = all.iter.lifetimes$lifetime[all.iter.lifetimes$deployment == j] - 1)
    zero.vector[1:length(zero.vector)] <- 0
    zero.vector[length(zero.vector) + 1] <- 1
    
    # pack into df
    focal.df <- data.frame(
      
      j = j,
      t = 1:length(zero.vector),
      event = zero.vector
      
    )
    
    km.data <- rbind(km.data, focal.df)
    
  }
  
  # truncate data at 100 weeks
  km.data.1 <- km.data %>% filter(t < 101)
  
  # apply KM function
  km.curve.data <- kap_mei(x = km.data.1,
                           response = "event")
  
  # add in iteration
  km.curve.data$iter = i
  
  # pack into df
  all.km.curves <- rbind(all.km.curves, km.curve.data)
  
  # print status (every 100 iterations)
  
  if (i %% 100 == 0) {
    
    print(paste0("Completed iteration ", i, " of ", nrow(draws))) 
    
  }
  
}

#_______________________________________________________________________________________________
# 8. Plot with empirical curve ----

# here we'll use 2500 iterations
some.km.curves <- all.km.curves %>% 
  
  group_by(iter) %>%
  
  slice(sample(1:max(all.km.curves$iter), 2500))

#_______________________________________________________________________________________________

# plot
ggplot() +
  
  theme_bw() +
  
  # simulated curves
  geom_line(data = some.km.curves,
            aes(x = t,
                y = S,
                group = iter),
            color = "aquamarine3",
            linewidth = 0.25,
            alpha = 0.05) +
  
  # empirical curve confidence envelope
  geom_ribbon(data = km.df.3,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.15) +
  
  # empirical curve
  geom_line(data = km.df.3,
            aes(x = t,
                y = S)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(4.25, 70))

# 350 x 350

