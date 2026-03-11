# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 03 - Censor-as-dead scenarios
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 12 Nov 2025 
# Date completed: 12 Nov 2025 
# Date last modified: 03 Feb 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in and prepare data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_final_cleaned.csv")

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
# 3. Logistic regression model ----
#_______________________________________________________________________________________________

glm.model <- glm(y.mort ~ Collar.type.1 + Transmitter.lifetime,
                 data = fates.1,
                 family = binomial)

summary(glm.model)

#_______________________________________________________________________________________________
# 4. Predict ----
#_______________________________________________________________________________________________

glm.pred <- predict(glm.model, fates.to.pred, type = "response", se.fit = T)

# data frame
glm.pred.1 <- data.frame(
  
  deployment.1 = fates.to.pred$deployment.1,
  mean = glm.pred$fit,
  l95 = glm.pred$fit - 1.96 * glm.pred$se.fit,
  u95 = glm.pred$fit + 1.96 * glm.pred$se.fit
  
) %>%
  
  # allow for conditional formatting
  # here we'll give an index for each of the two "censor-as-dead" scenarios
  mutate(
    
    scenario = case_when(
      
      # scenario 2, only those above 75% prob get flipped
      l95 > 0.75 ~ "scen2",
      
      # scenario 3, those above 50% get flipped
      l95 > 0.5 & l95 <= 0.75 ~ "scen3",
      
      # any others don't get flipped at all
      l95 <= 0.5 ~ "scen0"
    )
    
  ) %>%
  
  # arrange by lower interval
  arrange(l95) %>%
  
  mutate(unique.arrange = 1:n())

#_______________________________________________________________________________________________
# 6. Visualize ----
#_______________________________________________________________________________________________

ggplot(data = glm.pred.1) +
  
  theme_bw() +
  
  geom_hline(yintercept = 0.5,
             linetype = "dashed") +
  
  geom_hline(yintercept = 0.75,
             linetype = "dashed") +
  
  geom_errorbar(aes(x = unique.arrange,
                    ymin = l95,
                    ymax = u95,
                    color = scenario),
                width = 0,
                linewidth = 0.55) +
  
  geom_point(aes(x = unique.arrange,
                 y = mean,
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

fates.to.pred.2 <- glm.pred.1 %>%
  
  # binary indicator for which ones get flipped
  mutate(
    
    # scenario 2 (75%)
    flip.scen2 = ifelse(l95 > 0.75,
                        1,
                        0),
    
    # scenario 3 (50%)
    flip.scen3 = ifelse(l95 > 0.5,
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
    study.week,
    Sex.1,
    Collar.type.1,
    p.dm,
    p.o,
    p.js,
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

write.csv(fates.new.3, "Cleaned data/fates_forModel.csv", row.names = F)
