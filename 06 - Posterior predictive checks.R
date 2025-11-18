# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 06 - Posterior predictive checks
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 17 Nov 2025
# Date last modified: 18 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 0. Explanation ----
#_______________________________________________________________________________________________

# Bayesian p-value
# Here we'll apply a simplified procedure outlined in Jones et al. 2020
# For each week of our study, we'll use posterior draws to generate 
# predictions for our Poisson intensities, then calculate the frequency of 
# mortalities from simulated Poisson draws, then compare to observed frequencies

# we'll calculate the conditional (i.e., by week) difference between simulated and observed
# sum up all of these differences,
# do this for all iterations
# and look at how many iterations have more (conditional) simulated events than observed events
# should be close to 0.5

# Simulated lifetimes
# Here we'll simulate lifetimes from a subset of posterior draws, create Kaplan-Meier
# survivorship curves for each "population", then compare to the empirical K-M curve
# I'll just need to snap every individual to the same starting place

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- read.csv("Model outputs/model_1.csv")

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
        study.week,
        study.year.week
        
      ) %>%
      
      rename(week = study.year.week) %>%
      
      group_by(year, week) %>%
      
      slice(1),
    
    by = c("year", "week")
    
  ) %>%
  
  # standardize BCI.1
  mutate(BCI.s = (BCI.1 - mean(BCI.1)) / sd(BCI.1))

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
# 4. Split dataset by study.week ----
#_______________________________________________________________________________________________

fates.1.split <- split(fates.1, fates.1$study.week)

#_______________________________________________________________________________________________
# 5. Write Bayesian p-value function ----

# lapply() was a good thought, but I think just looping by iteration and doing the summation
# within the function is probably fine. We can also print the status by iteration

# helper function 
# functions
calc_intensity <- function (y) {
  
  # BLH spline prediction
  # weights w
  w <- t(as.matrix(as.numeric(c(iter.draw[paste0("w", ".", y$sex_forest, "..1.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..2.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..3.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..4.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..5.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..6.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..7.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..8.")],
                                iter.draw[paste0("w", ".", y$sex_forest, "..9.")]))))
  
  # multiply and sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
  
  # add to intercept and exponentiate
  blh = as.numeric(exp(iter.draw[paste0("a0.", y$sex_forest, ".")] + w.by.b.sum[y$week]))
  
  # coefficient prediction (must be on the log scale!)
  coef.pred = exp(log(iter.draw$hr_bci) * y$BCI.s +
                    log(iter.draw$hr_pil_total1) * y$post1 * y$pil +
                    log(iter.draw$hr_pil_total2) * y$post2 * y$pil +
                    log(iter.draw$hr_ret_total1) * y$post1 * y$ret +
                    log(iter.draw$hr_ret_total2) * y$post2 * y$ret)
  
  # total intensity
  intens = blh * coef.pred
  
  # Poisson draw
  pois.draw = rpois(n = 1, intens)
  
  # pack into a df
  y.intens.pois <- data.frame(intens = intens,
                              pois.draw = ifelse(pois.draw > 0, 1, pois.draw))
  
  return(y.intens.pois)
  
}

#_______________________________________________________________________________________________

bayes_p_mort <- function(
    
  draws = model.fit.1,
  response = "y.mort.scen1"
  
  ) {
  
  # loop through iterations
  discrep.sim.all.draws <- vector(length = nrow(draws))
  discrep.obs.all.draws <- vector(length = nrow(draws))
  
  for (i in 1:nrow(draws)) {
    
    # iteration
    iter.draw <- draws[i, ]
    
    # loop through study.weeks
    discrep.sim <- vector(length = max(fates.1$study.week))
    discrep.obs <- vector(length = max(fates.1$study.week))
    
    for (j in 1:max(fates.1$study.week)) {
      
      focal.week <- fates.1[fates.1$study.week == j, ]
      
      # function to apply through all observations
      indiv.split <- split(focal.week, focal.week$deployment.1)
      
      # apply function
      sim.events <- do.call(rbind, lapply(indiv.split, calc_intensity))
      
      # determine the expected number of events
      # this is just the sum of all intensities
      e.events = sum(sim.events$intens)
      
      # calculate the conditional discrepancy
      # Dj = ((sim/obs events - predicted events)^2) / predicted events
      discrep.sim[j] = ((sum(sim.events$pois.draw) - sum(sim.events$intens))^2) / sum(sim.events$intens)
      discrep.obs[j] = ((sum(focal.week[[response]]) - sum(sim.events$intens))^2) / sum(sim.events$intens)
      
      
      
    }
    
    # sum over all study.weeks
    discrep.sim.all.draws[i] = sum(discrep.sim)
    discrep.obs.all.draws[i] = sum(discrep.obs)
      
    # print status (every 100 iterations)
    
    if (i %% 5 == 0) {
      
      print(paste0("Completed iteration ", i, " of ", nrow(draws))) 
      
    }
    
  }
  
  # bind together
  discrep.sim.obs <- data.frame(sim = discrep.sim.all.draws,
                                obs = discrep.obs.all.draws)
  
  return(discrep.sim.obs)

}

#_______________________________________________________________________________________________
# 6. Apply function ----

# for all posterior iterations, this would take a long time. Best to leave it overnight

#_______________________________________________________________________________________________

bayes.diff <- bayes_p_mort(draws = draws, response = "y.mort.scen1")


bayes.diff <- data.frame(
  
  sim = discrep.sim.all.draws,
  obs = discrep.obs.all.draws,
  diff = discrep.sim.all.draws - discrep.obs.all.draws
  
  )

#_______________________________________________________________________________________________
# 7. Plot and calculate Bayes p-value ----
#_______________________________________________________________________________________________

length(which(bayes.diff$diff > 0)) / nrow(bayes.diff) # 0.12

# distribution
ggplot(data = bayes.diff,
       aes(x = diff)) +
  
  theme_bw() +
  
  geom_density(fill = "gray",
               alpha = 0.25) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  
  xlab("Dsim - Dobs")

# 1:1 plot
ggplot(data = bayes.diff,
       aes(x = obs,
           y = sim)) +
  
  theme_bw() +
  
  geom_point(shape = 21) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  
  xlab("Observed discrepancy") +
  ylab("Simulated discrepancy") +
  
  coord_cartesian(xlim = c(100, 250),
                  ylim = c(100, 250))

#_______________________________________________________________________________________________
# 8. Prepare dataset for Kaplan-Meier curves ----
#_______________________________________________________________________________________________
# 8a. Prepare empirical data ----

# Here we'll keep everything the same, but include a "t" variable denoting the week from the beginning of 
# monitoring per deployment

#_______________________________________________________________________________________________

# split by deployment
fates.1.deploy.split <- split(fates.1, fates.1$deployment.1)

# function to incorporate t since monitoring started
add_t <- function (x) {
  
  x$t = 1:nrow(x)
  
  return(x)
  
}

# apply
fates.2 <- do.call(rbind, lapply(fates.1.deploy.split, add_t))

#_______________________________________________________________________________________________
# 8b. Write function ----

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

km.df <- kap_mei(fates.2)

# plot
ggplot(data = km.df,
       aes(x = t,
           y = S)) +
  
  theme_bw() +
  
  geom_ribbon(aes(ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.25) +
  
  geom_line() +
  
  theme(panel.grid = element_blank()) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1))

#_______________________________________________________________________________________________
# 9. Simulate lifetimes ----

# Here we'll loop through each iteration (take a subsample first),
# loop through each deployment (n = 383),
# applying a calc intensity function on an 86-week dataset,
# take a Poisson draw and truncate the dataset after the first event
# then apply the K-M curve function

# the output will be a collection of K-M curves

# new calc_intensity function
calc_intensity_2 <- function (y) {
  
  # BLH spline prediction
  # weights w
  w <- t(as.matrix(as.numeric(c(iter.draw[paste0("w", ".", y$sex_forest[1], "..1.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..2.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..3.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..4.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..5.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..6.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..7.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..8.")],
                                iter.draw[paste0("w", ".", y$sex_forest[1], "..9.")]))))
  
  # multiply and sum to create "normal" scale spline prediction
  w.by.b.sum <- apply(sweep(basis.pred, 2, w, `*`), 1, sum)
  
  # add to intercept and exponentiate
  blh = as.numeric(exp(iter.draw[paste0("a0.", y$sex_forest, ".")] + w.by.b.sum[y$week]))
  
  # coefficient prediction (must be on the log scale!)
  coef.pred = exp(log(iter.draw$hr_bci) * y$BCI.s +
                    log(iter.draw$hr_pil_total1) * y$post1 * y$pil +
                    log(iter.draw$hr_pil_total2) * y$post2 * y$pil +
                    log(iter.draw$hr_ret_total1) * y$post1 * y$ret +
                    log(iter.draw$hr_ret_total2) * y$post2 * y$ret)
  
  # total intensity
  intens = blh * coef.pred
  
  # Poisson draw
  pois.draw = rpois(n = 1, intens)
  
  # pack into a df
  y.intens.pois <- data.frame(intens = intens,
                              pois.draw = rpois(length(intens), intens))
  
  return(y.intens.pois)
  
}

#_______________________________________________________________________________________________

sim_lifetime <- function (
    
  x = fates.2,
  draws = model.fit.1
  
) {
  
  # loop through draws
  all.km.curves <- data.frame()
  
  for (i in 1:nrow(draws)) {
    
    iter.draw = draws[i, ]
    
    # apply to each deployment
    deploy.split <- split(x, x$deployment.1)
    
    # function
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
      deploy.df <- cbind(deploy.df, calc_intensity_2(y = deploy.df))
      
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
    
    # apply function
    all.iter.lifetimes <- do.call(rbind, lapply(deploy.split, sim_lifetime_byDeploy))
    
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
    
    if (i %% 5 == 0) {
      
      print(paste0("Completed iteration ", i, " of ", nrow(draws))) 
      
    }

  }
  
  return(all.km.curves)
  
}

#_______________________________________________________________________________________________
# 9a. Apply function ----
#_______________________________________________________________________________________________

km.curves <- sim_lifetime(x = fates.2,
                          draws = model.fit.1[sample(1:15000, size = 50), ])

#_______________________________________________________________________________________________
# 10. Plot with empirical curve ----
#_______________________________________________________________________________________________

# plot
ggplot() +
  
  theme_bw() +
  
  # simulated curves
  geom_line(data = km.curves,
            aes(x = t,
                y = S,
                group = iter),
            color = "aquamarine3",
            linewidth = 0.25) +
  
  # empirical curve confidence envelope
  geom_ribbon(data = km.df,
              aes(x = t,
                  y = S,
                  ymin = ci.low,
                  ymax = ci.upp),
              alpha = 0.15) +
  
  # empirical curve
  geom_line(data = km.df,
            aes(x = t,
                y = S)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  
  xlab("Weeks since deployment") +
  
  ylab("Cumulative survival") +
  
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(4.25, 70))
