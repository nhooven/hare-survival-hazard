# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Censor-as-dead scenarios
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Nov 2025 
# Date completed: 12 Nov 2025 
# Date last modified: 14 Nov 2025 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(nimble)          # MCMC
library(coda)            # play nice with posterior samples
library(bayestestR)

#_______________________________________________________________________________________________
# 2. Read in and prepare data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_final_cleaned_2.csv")

# subset to those observations we'll fit the model to
fates.1 <- fates %>% 
  
  filter((Event.type == "Mortality" |
         (Event.type == "Censor" & General.cause == "Dead transmitter"))) %>%
  
  group_by(deployment.1) %>%
  
  slice(n()) %>%
  
  ungroup() %>%
  
  # mort variable
  mutate(y.mort = ifelse(Event.type == "Mortality",
                         1,
                         0)) %>%
  
  dplyr::select(Collar.type.1,
                Transmitter.lifetime,
                y.mort)

# dataset to predict on
fates.to.pred <- fates %>%
  
  filter(Event.type == "Censor" &
         General.cause == "Unknown" &
         event == 1)

#_______________________________________________________________________________________________
# 3. Prepare for modeling ----
#_______________________________________________________________________________________________

constants = list(
  
  N = nrow(fates.1),                   # observations in fates df for modeling
  N_to_pred = nrow(fates.to.pred)      # observations to predict upon
  
)

fates.data <- list(
  
  # for model fitting
  collar = fates.1$Collar.type.1,
  tag_life = fates.1$Transmitter.lifetime,
  y_mort = fates.1$y.mort,
  
  # for prediction
  collar_to_pred = fates.to.pred$Collar.type.1,
  tag_life_to_pred = fates.to.pred$Transmitter.lifetime
  
)

#_______________________________________________________________________________________________
# 4. Set up model ----
#_______________________________________________________________________________________________
# 4a. Code ----
#_______________________________________________________________________________________________

model.code <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0, sd = 1)
  b1 ~  dt(0, sigma = 2.5, df = 1)        # collar
  b2 ~  dt(0, sigma = 2.5, df = 1)        # tag life
  
  # likelihood
  for (i in 1:N) {
    
    y_mort[i] ~ dbern(p_mort[i])
    
    logit(p_mort[i]) <- b0 + b1 * collar[i] + b2 * tag_life[i]
    
  }
  
  # prediction
  for (k in 1:N_to_pred) {
    
    p_pred[k] <- ilogit(b0 + b1 * collar_to_pred[k] + b2 * tag_life_to_pred[k])
    
  }
  
})

# parameters monitored
params <- c("b0",
            "b1",
            "b2",
            "p_pred")

# inits
inits <- list(
  
  b0 = rnorm(1, 0, 2.5),
  b1 = rnorm(1, 0, 2.5),
  b2 = rnorm(1, 0, 2.5)
  
)

#_______________________________________________________________________________________________
# 4b. Model ----
#_______________________________________________________________________________________________

model.pred <- nimbleModel(code = model.code,
                          constants = constants,
                          data = fates.data,
                          inits = inits) 

#_______________________________________________________________________________________________
# 4c. Run MCMC ----
#_______________________________________________________________________________________________

# run model
model.fit <- nimbleMCMC(model = model.pred,
                        monitors = params,
                        nchains = 3,
                        nburnin = 5000,
                        niter = 20000,
                        thin = 1,
                        samplesAsCodaMCMC = TRUE)

summary(model.fit)

HPDinterval(model.fit)

effectiveSize(model.fit)

#_______________________________________________________________________________________________
# 5. Extract HDIs for visualization ----
#_______________________________________________________________________________________________

# pack all chains together
model.samples <- do.call(rbind, model.fit)

# keep only the predictions as a df
model.samples.1 <- as.data.frame(model.samples[ , 4:ncol(model.samples)])

# calculate HDIs
model.hdi <- hdi(model.samples.1, ci = 0.95)

# calculate posterior median
model.median <- apply(model.samples.1, 2, median)

# put back into data.frame
fates.to.pred.1 <- fates.to.pred %>%
  
  mutate(med = model.median,
         low.hdi = model.hdi$CI_low,
         upp.hdi = model.hdi$CI_high) %>%
  
  # allow for conditional formatting
  # here we'll give an index for each of the two "censor-as-dead" scenarios
  mutate(
    
    scenario = case_when(
      
      # scenario 2, only those above 75% prob get flipped
      low.hdi > 0.75 ~ "scen2",
      
      # scenario 3, those above 50% get flipped
      low.hdi > 0.5 & low.hdi <= 0.75 ~ "scen3",
      
      # any others don't get flipped at all
      low.hdi <= 0.5 ~ "scen0"
    )
    
  ) %>%
  
  # arrange by lower interval
  arrange(low.hdi) %>%
  
  mutate(unique.arrange = 1:n())

#_______________________________________________________________________________________________
# 6. Visualize ----
#_______________________________________________________________________________________________

ggplot(data = fates.to.pred.1) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0.5,
             linetype = "dashed") +
  
  geom_hline(yintercept = 0.75,
             linetype = "dashed") +
  
  geom_errorbar(aes(x = unique.arrange,
                    ymin = low.hdi,
                    ymax = upp.hdi,
                    color = scenario),
                width = 0,
                linewidth = 0.55) +
  
  geom_point(aes(x = unique.arrange,
                 y = med,
                 fill = scenario),
             shape = 21,
             size = 0.55) +
  
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = c(0.6, 0.2),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  
  ylab("Probability of mortality event") +
  
  scale_fill_manual(values = c("darkgray", "purple", "orange"),
                    labels = c("Always censored",
                               "Mortality in scenarios 2 and 3",
                               "Mortality in scenario 3")) +
  scale_color_manual(values = c("darkgray", "purple", "orange"),
                     labels = c("Always censored",
                                "Mortality in scenarios 2 and 3",
                                "Mortality in scenario 3"))

# 380 x 350

#_______________________________________________________________________________________________
# 7. Flip outcomes in new dataset ----
#_______________________________________________________________________________________________

fates.to.pred.2 <- fates.to.pred.1 %>%
  
  # binary indicator for which ones get flipped
  mutate(
    
    # scenario 2 (75%)
    flip.scen2 = ifelse(low.hdi > 0.75,
                        1,
                        0),
    
    # scenario 3 (50%)
    flip.scen3 = ifelse(low.hdi > 0.5,
                        1, 
                        0)
    
  ) %>%
  
  # keep only columns we need to create new vars in full dataset
  dplyr::select(deployment.1,
                flip.scen2,
                flip.scen3)

# clone y.mort / y.pred columns in original dataset
fates.new <- fates

fates.new <- fates.new %>%
  
  mutate(y.mort.scen2 = y.mort.scen1,
         y.mort.scen3 = y.mort.scen1,
         y.pred.scen2 = y.pred.scen1,
         y.pred.scen3 = y.pred.scen1)

# write function that will "flip" zeroes to ones according to fates.to.pred.2
flip_fates <- function (df) {

  if (df$deployment.1[1] %in% fates.to.pred.2$deployment.1) {
    
    # index of event row
    event.ind = which(df$event == 1)
    
    # flip if needed
    df$y.mort.scen2[event.ind] = fates.to.pred.2$flip.scen2[fates.to.pred.2$deployment.1 == df$deployment.1[1]]
    df$y.pred.scen2[event.ind] = fates.to.pred.2$flip.scen2[fates.to.pred.2$deployment.1 == df$deployment.1[1]]
    
    df$y.mort.scen3[event.ind] = fates.to.pred.2$flip.scen3[fates.to.pred.2$deployment.1 == df$deployment.1[1]]
    df$y.pred.scen3[event.ind] = fates.to.pred.2$flip.scen3[fates.to.pred.2$deployment.1 == df$deployment.1[1]]

  }
  
  return(df)
  
}

# split by deployment
fates.new.split <- split(fates.new, f = fates.new$deployment.1)

# apply function
fates.new.2 <- do.call(rbind, lapply(fates.new.split, flip_fates))

# I think this worked, let's check
sum(fates.new.2$y.mort.scen1)
sum(fates.new.2$y.mort.scen2)
sum(fates.new.2$y.mort.scen3)

sum(fates.new.2$y.pred.scen1)
sum(fates.new.2$y.pred.scen2)
sum(fates.new.2$y.pred.scen3)

# excellent, now let's clean this up and get ready for modeling

#_______________________________________________________________________________________________
# 8. Clean up ----
#_______________________________________________________________________________________________

fates.new.3 <- fates.new.2 %>%
  
  # only columns we need, in the right order
  dplyr::select(
    
    deployment.1,
    cluster,
    site,
    sex_forest,
    MRID,
    Event.type,
    General.cause,
    y.mort.scen1,
    y.mort.scen2,
    y.mort.scen3,
    y.pred.scen1,
    y.pred.scen2,
    y.pred.scen3,
    week,
    Sex.1,
    Collar.type.1,
    BCI.1,
    BCI.2,
    post1,
    post2,
    ret,
    pil
    
  ) %>%
  
  # rename as needed
  rename(sex = Sex.1,
         collar = Collar.type.1) %>%
  
  # sex_forest as integer
  mutate(sex_forest = as.integer(as.factor(sex_forest)))

#_______________________________________________________________________________________________
# 9. Write to .csv ----
#_______________________________________________________________________________________________

write.csv(fates.new.3, "Cleaned data/fates_forModel.csv")
