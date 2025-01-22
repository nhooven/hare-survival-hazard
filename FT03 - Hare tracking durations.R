# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Fate timing and incidence
# Script: FT03 - Hare tracking durations
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Jan 2025
# Date completed: 22 Jan 2025
# Date last modified: 21 Jan 2025
# R version: 4.2.2

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(lubridate)       # work with dates

#_______________________________________________________________________________________________
# 2. Read in data ----
#_______________________________________________________________________________________________

fates <- read.csv("Cleaned data/fates_timing_cleaned_01_21_2025.csv")

# day lookup table (for reference)
day.lookup <- read.csv("Cleaned data/day_lookup.csv")

#_______________________________________________________________________________________________
# 3. Monthly cutoffs ----
#_______________________________________________________________________________________________

month.cutoffs <- data.frame(breaks = c(1, 5.57, 9.86, 
                                       14.29, 18.71, 22.71, 
                                       27.14, 31.43, 35.86, 
                                       40.14, 44.57, 49.86),
                            labels = c("Oct", "Nov", "Dec", 
                                       "Jan", "Feb", "Mar", 
                                       "Apr", "May", "Jun", 
                                       "Jul", "Aug", "Sep"))

#_______________________________________________________________________________________________
# 4. Create and clean variables ----
#_______________________________________________________________________________________________
# 4a. Generic "event" variable ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  mutate(Event = case_when(Event.type == "Censor" & General.cause == "Removed" ~ "Removed",
                           Event.type == "Censor" & General.cause != "Removed" ~ "Censor",
                           Event.type == "On air" ~ "",
                           Event.type == "Mortality" ~ "Mortality")) %>%
  
  # and reorder
  mutate(Event = factor(Event,
                        levels = c("Censor", "Mortality", "Removed", "")))

#_______________________________________________________________________________________________
# 4b. Date format ----
#_______________________________________________________________________________________________

fates.1 <- fates.1 %>% 
  
  mutate(Estimated.event.date = as.Date(Estimated.event.date),
         Capture.date = as.Date(Capture.date))

#_______________________________________________________________________________________________
# 5. Plots ----
#_______________________________________________________________________________________________
# 5a. All - why not? ----
#_______________________________________________________________________________________________

# order by Capture.date (within each treatment type)
fates.2 <- fates.1 %>%
  
  arrange(Capture.date) %>%
  
  group_by(trt) %>%
  
  mutate(identifier = 1:n(),
         
         Cluster = factor(cluster, 
                          levels = c("1", "2", "3", "4")))

ggplot() +
  
  theme_bw() +
  
  # add treatment period as a gray rectangle
  annotate("rect", 
           xmin = as.Date("2023-10-04"),
           xmax = as.Date("2023-10-12"),
           ymin = Inf,
           ymax = -Inf,
           fill = "lightgray") +
  
  facet_wrap(~ trt,
             ncol = 1,
             scales = "free_y") +
  
  geom_segment(data= fates.2,
               aes(x = Capture.date,
                   xend = Estimated.event.date,
                   y = identifier,
                   group = Animal.ID,
                   color = Cluster)) +
  
  # segment color
  scale_color_viridis_d(begin = 0,
                        end = 0.8) +
  
  geom_point(data = fates.2,
             aes(y = identifier,
                 group = Animal.ID,
                 x = Estimated.event.date,
                 shape = Event,
                 fill = Event,
                 size = Event)) +
  
  # point fills
  scale_fill_manual(values = c("yellow", "red", "white", NA)) +
  
  # point shapes
  scale_shape_manual(values = c(21, 22, 23, NA)) +
  
  # point sizes 
  scale_size_manual(values = c(1.0, 0.9, 0.9, NA)) +
  
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(hjust = 0)) +
  
  xlab("Date") +
  
  scale_y_reverse()
