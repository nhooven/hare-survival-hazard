# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 07 - Schoenfeld residuals
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 14 Nov 2025 
# Date completed: 17 Nov 2025 
# Date last modified: 05 Nov 2026 
# R version: 4.5.2

#_______________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________

# model samples
model.fit.1 <- readRDS("models/model_1.rds")

# dataset
fates <- read.csv("data_cleaned/fates_forModel.csv")

# previous dataset because it has the year variable
fates.year <- read.csv("data_cleaned/fates_final_cleaned_2.csv")

# BCI by deployment
fates.deploy <- readRDS("data_cleaned/fates_deploy.rds")
fates.deploy.df <- fates.deploy$fates.deploy

# day lookup table
day.lookup <- read.csv("data_cleaned/day_lookup.csv")

#_______________________________________________________________________________
# 3. Prepare data ----
#_______________________________________________________________________________
# 3a. Betas ----

# function to turn mcmc.list into df and select only columns we need
extract_b <- function (x) {
  
  x.1 <- do.call(rbind, x)
  
  x.1.hr <- as.data.frame(x.1) %>%
    
    dplyr::select(
      
      b_bci,
      b_col,
      b_dm,
      b_open,
      b_post1,
      b_post2,
      b_ret_post1,
      b_ret_post2,
      b_pil_post1,
      b_pil_post2
      
    )
  
  return(x.1.hr)
  
}

#_______________________________________________________________________________

model.fit.1.b <- extract_b(model.fit.1)

#_______________________________________________________________________________
# 3b. Attribute study-week to fates df ----
#_______________________________________________________________________________

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
  
  # join in BCI data
  left_join(fates.deploy.df %>% dplyr::select(deployment.1, BCI)) %>%
  
  # impute mean BCI 
  replace_na(list(BCI = fates.deploy$bci.mean)) %>%
  
  # standardize variables
  mutate(BCI.s = (BCI - fates.deploy$bci.mean) / fates.deploy$bci.sd,
         p.dm.s = (p.dm - mean(p.dm)) / sd(p.dm),
         p.open.s = (p.o - mean(p.o)) / sd(p.o))

#_______________________________________________________________________________
# 3c. Split dataset by study year ----
#_______________________________________________________________________________

fates.1.pre <- fates.1 %>% filter(post1 == 0 & post2 == 0)
fates.1.post1 <- fates.1 %>% filter(post1 == 1 & post2 == 0)
fates.1.post2 <- fates.1 %>% filter(post1 == 0 & post2 == 1)

#_______________________________________________________________________________
# 4. Function - calculate scaled Schoenfeld residuals ----

# we'll do this for BCI, COL, TRT, and landscape variables

#_______________________________________________________________________________

calc_schoen_resid <- function (
    
  .coefs,
  .fates
  
) {
  
  # calculate VCOV matrix
  vcov <- matrix(data = NA, nrow = ncol(.coefs), ncol = ncol(.coefs))
  
  # sequence
  coef.seq <- 1:ncol(.coefs)
  
  # loop through
  for (i in coef.seq) {
    
    # variance
    vcov[i, i] <- var(.coefs[ , i])
    
    # covariances
    for (j in coef.seq[coef.seq != i]) {
      
      vcov[i, j] <- cov(.coefs[ , i], .coefs[ , j])
      vcov[j, i] <- cov(.coefs[ , i], .coefs[ , j])
      
    } # j
    
  } # i
  
  # inverse of the VCOV
  vcov.inv <- solve(vcov)
  
  # names
  rownames(vcov) <- colnames(.coefs)
  colnames(vcov) <- colnames(.coefs)
  
  # only proceed if there are > 0 events
  if (1 %in% .fates[["y.mort1"]]) {
   
    # extract covariate values for all individuals at risk
    x.vals <- .fates |> dplyr::select(
      
      BCI.s, 
      collar,
      p.dm.s,
      p.open.s,
      post1,
      post2,
      ret,
      pil
      
      ) |>
      
      # "interaction" covariates
      mutate(ret.post1 = ret * post1,
             ret.post2 = ret * post2,
             pil.post1 = pil * post1,
             pil.post2 = pil * post2) |>
      
      # keep columns
      dplyr::select(
        
        BCI.s, 
        collar,
        p.dm.s,
        p.open.s,
        post1,
        post2,
        ret.post1,
        ret.post2,
        pil.post1,
        pil.post2
        
      )
    
    # calculate risk-weighted average for all covariates
    # sum both numerator and denominator over all individuals at risk
    # loop over posterior samples
    rw.num.all <- array(NA, dim = c(nrow(.fates),
                                    ncol(x.vals),
                                    nrow(.coefs)))
    
    rw.denom.all <- array(NA, dim = c(nrow(.fates),
                                      ncol(x.vals),
                                      nrow(.coefs)))
    
    for (i in 1:nrow(.coefs)) {
      
      # define matrix and vector
      x.mat <- as.matrix(x.vals)
      coef.vec <- as.numeric(.coefs[i, ])
      
      # num/denom for summing
      rw.num.all[ , , i] <- x.mat * exp(sweep(x.mat, 2, coef.vec, `*`))
      rw.denom.all[ , , i] <- exp(sweep(x.mat, 2, coef.vec, `*`))
      
    } # i - n posterior samples
    
    # sum (rowwise) and divide
    rw.num.sum <- apply(rw.num.all, c(3, 2), sum)
    rw.denom.sum <- apply(rw.denom.all, c(3, 2), sum)
    
    # yielding a n.samples x n.covariates matrix
    rw <- rw.num.sum / rw.denom.sum
    
    # calculate the residual(s) for individual(s) experiencing the event
    x.event <- .fates |>
      
      filter(y.mort1 == 1) |>
      
      dplyr::select(
        
        BCI.s, 
        collar,
        p.dm.s,
        p.open.s,
        post1,
        post2,
        ret,
        pil
        
      ) |>
      
      # "interaction" covariates
      mutate(ret.post1 = ret * post1,
             ret.post2 = ret * post2,
             pil.post1 = pil * post1,
             pil.post2 = pil * post2) |>
      
      # keep columns
      dplyr::select(
        
        BCI.s, 
        collar,
        p.dm.s,
        p.open.s,
        post1,
        post2,
        ret.post1,
        ret.post2,
        pil.post1,
        pil.post2
        
      )
    
    # loop through events to calculate resid vector for each covariate
    r.ik.all <- data.frame()
    
    for (i in 1:nrow(x.event)) {
      
      x.row <- as.numeric(x.event[i, ])
      
      # posterior of Schoenfeld residuals
      r.ik <- sweep(rw, 2, x.row, `-`)
      
      # posterior of scaled Schoenfeld residuals
      # assume that the VCOV is constant in time
      # https://stats.stackexchange.com/questions/547078/schoenfeld-residuals-plain-english-explanation-please
      r.star.ik <- matrix(NA, nrow = nrow(.coefs), ncol = ncol(r.ik))
      
      for (j in 1:nrow(.coefs)) {
        
        coef.vec <- as.numeric(.coefs[j, ])
        r.vec <- as.numeric(r.ik[j, ])
        
        r.star.ik[j, ] <- coef.vec + vcov.inv %*% r.vec
        
      } # j - posterior samples
      
      # means, SD, 90% CI
      r.ik.sum <- data.frame(
        
        indiv.week = paste0(i, "_", .fates$study.week[1]),
        study.week = .fates$study.week[1],
        cov = colnames(.coefs),
        mean = apply(r.star.ik, 2, mean),
        sd = apply(r.star.ik, 2, sd),
        l90 = apply(r.star.ik, 2, quantile, prob = 0.05),
        u90 = apply(r.star.ik, 2, quantile, prob = 0.95)
        
      )
      
      r.ik.all <- rbind(r.ik.all, r.ik.sum)
      
    } # i - events == 1
     
    return(r.ik.all)
    
  }
  
} # fn

#_______________________________________________________________________________
# 4b. Use function ----
#_______________________________________________________________________________

schoen.resid <- do.call(rbind, 
                        lapply(split(fates.1, fates.1$study.week), 
                               calc_schoen_resid,
                               .coefs = model.fit.1.b))

#_______________________________________________________________________________
# 5. Plots ----
#_______________________________________________________________________________
# 5a. Prepare for plotting ----
#_______________________________________________________________________________

# add post-treatment periods
first.post1 <- min(fates.1$study.week[fates.1$post1 == 1])
first.post2 <- min(fates.1$study.week[fates.1$post2 == 1])

schoen.resid <- schoen.resid |>
  
  mutate(period = case_when(
    
    study.week < first.post1 ~ "PRE",
    study.week >= first.post1 & study.week < first.post2 ~ "POST1",
    study.week >= first.post2 ~ "POST2"
    
  )
  
  )

# split datasets
# continuous covariates
schoen.resid.continuous <- schoen.resid |> 
  
  filter(cov %in% c("b_bci",
                    "b_col",
                    "b_dm",
                    "b_open")) |>
  
  # factor levels
  mutate(cov = factor(cov,
                      labels = c("BCI", "COLLAR", "DM", "OPEN")))

# treatment covariates
schoen.resid.trt <- schoen.resid |>
  
  filter(period %in% c("POST1", "POST2")) |>
  
  filter(cov %in% c("b_post1",
                    "b_post2",
                    "b_ret_post1",
                    "b_ret_post2",
                    "b_pil_post1",
                    "b_pil_post2"))

# split and recombine
# we need to do this because the residuals for YR1 in YR2 are nonsensical
schoen.resid.trt.1 <- bind_rows(
  
  # YR 1
  schoen.resid.trt |>
    
    filter(cov %in% c("b_post1",
                      "b_ret_post1",
                      "b_pil_post1") &
           period == "POST1"),
  
  # YR 2
  schoen.resid.trt |>
    
    filter(cov %in% c("b_post2",
                      "b_ret_post2",
                      "b_pil_post2") &
           period == "POST2")
  
) |>
  
  # new covariate level
  mutate(cov1 = case_when(
    
    cov %in% c("b_post1", "b_post2") ~ "YR",
    cov %in% c("b_ret_post1", "b_ret_post2") ~ "RET",
    cov %in% c("b_pil_post1", "b_pil_post2") ~ "PIL"
    
  )) |>
  
  # factor order
  mutate(cov1 = factor(cov1, levels = c("YR", "RET", "PIL")))

#_______________________________________________________________________________
# 5b. Continuous covariates ----
#_______________________________________________________________________________

ggplot(data = schoen.resid.continuous) +
  
  theme_bw() +
  
  facet_wrap(~ cov,
             scales = "free_y") +
  
  geom_errorbar(aes(x = study.week,
                    y = mean,
                    ymin = l90,
                    ymax = u90),
                width = 0) +
  
  geom_point(aes(x = study.week,
                 y = mean),
             size = 0.5) +
  
  geom_smooth(aes(x = study.week,
                  y = mean,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA),
        legend.position = "none") +
  
  xlab("Study week") +
  ylab("Scaled Schoenfeld residual")

# 500 x 500

#_______________________________________________________________________________
# 5c. Treatment covariates ----
#_______________________________________________________________________________

ggplot(data = schoen.resid.trt.1) +
  
  theme_bw() +
  
  facet_wrap(~ cov1,
             scales = "free_y") +
  
  geom_errorbar(aes(x = study.week,
                    y = mean,
                    ymin = l90,
                    ymax = u90),
                width = 0) +
  
  geom_point(aes(x = study.week,
                 y = mean),
             size = 0.5) +
  
  # POST1
  geom_smooth(data = schoen.resid.trt.1 |> filter(period == "POST1"),
              aes(x = study.week,
                  y = mean,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  # POST1
  geom_smooth(data = schoen.resid.trt.1 |> filter(period == "POST2"),
              aes(x = study.week,
                  y = mean,
                  weight = 1 / sd),
              method = "gam",
              se = T,
              color = "aquamarine4",
              fill = "aquamarine4") +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA),
        legend.position = "none") +
  
  xlab("Study week") +
  ylab("Scaled Schoenfeld residual")

# 750 x 286