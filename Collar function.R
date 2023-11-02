# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: Collar function
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 02 Nov 2023
# Date completed: 
# Date last modified: 
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

transmitters <- read.csv("transmitter_timelines.csv")

#_______________________________________________________________________________________________
# 3. Clean data for summaries and analysis ----
#_______________________________________________________________________________________________

transmitters.1 <- transmitters %>% 
  
  # remove all on-air transmitters
  filter(Status != "On air" &
           
  # and all transmitters that were never deployed         
         Total.days != 0 &
    
  # as well as transmitters whose monitoring was missed  
         Notes != "Left off sheet- do not include in calculation")

# subset all lost transmitters
transmitters.lost <- transmitters.1 %>%
  
  filter(Status == "Lost")

#_______________________________________________________________________________________________
# 4. Lost collar lifetimes ----
#_______________________________________________________________________________________________

# group by type
transmitters.lost %>% 
  group_by(Collar.type) %>%
  
  # summarize lifetimes
  summarize(n = n(),
            min = min(Total.days),
            max = max(Total.days),
            mean = mean(Total.days),
            median = median(Total.days))

# plot densities
ggplot(data = transmitters.lost,
       
       # separate colors for collar type
       aes(x = Total.days,
           color = Collar.type,
           fill = Collar.type)) +
  
  # theme
  theme_bw() +
  
  # density plots
  geom_density(alpha = 0.25,
               linewidth = 1.5)

#_______________________________________________________________________________________________
# 5. Probability of becoming lost ----
#_______________________________________________________________________________________________

transmitters.2 <- transmitters.1 %>% 
  
  # add binary variable for "lost/not lost"
  mutate(lost = ifelse(Status == "Lost",
                       1,
                       0))

# plot
ggplot(data = transmitters.2,
       aes(x = Total.days,
           y = lost,
           group = Collar.type,
           color = Collar.type,
           fill = Collar.type)) +
  
  # theme
  theme_bw() +
  
  # facet
  facet_wrap(~ Collar.type,
             scale = "free") +
  
  # add points
  geom_point() +
  
  # add logistic curve
  stat_smooth(method = "glm",  
              method.args = list(family = binomial)) +
  
  # remove legend
  theme(legend.position = "none")
