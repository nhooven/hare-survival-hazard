# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Fate timing and incidence
# Script: FT02 - Fate incidence bar charts
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 21 Jan 2025
# Date completed: 21 Jan 2025
# Date last modified: 22 Jan 2025
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

# remove "on air" collars
fates <- fates %>% filter(Event.type != "On air")

#_______________________________________________________________________________________________
# 3. Split data ----
#_______________________________________________________________________________________________

# censoring
fates.censor <- fates %>% filter(Event.type == "Censor" &
                                 General.cause != "Removed")

# predation mortality
fates.pred <- fates %>% filter(Event.type == "Mortality" & 
                               General.cause == "Predation")

# other mortality
fates.oth <- fates %>% filter(Event.type == "Mortality" & 
                              General.cause != "Predation")

#_______________________________________________________________________________________________
# 4. Weekly bar charts ----
#_______________________________________________________________________________________________
# 4a. Monthly cutoffs ----
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
# 4b. Censoring ----
#_______________________________________________________________________________________________

# reorder and label factor levels
fates.censor$Cause <- factor(fates.censor$General.cause,
                             levels = c("Dead transmitter",
                                        "Unknown",
                                        "Slipped",
                                        "Capture mortality"),
                             labels = c("Dead transmitter",
                                        "Unknown",
                                        "Slipped",
                                        "Collar-related mortality"))

ggplot(data = fates.censor,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1) +
  
  # theme arguments
  theme(
        # remove gridlines
        panel.grid = element_blank(),
        
        # adjust axis text
        axis.text.x = element_text(angle = 270,
                                   vjust = 0.25),
        
        # remove x axis title
        axis.title.x = element_blank(),
        
        # legend position
        legend.position = "top",
        
        # legend size
        legend.key.size = unit(0.35, "cm")
        ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.24, 5)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:5)

#_______________________________________________________________________________________________
# 4c. Predation mortality ----
#_______________________________________________________________________________________________

# reclassify
fates.pred <- fates.pred %>%
  
  mutate(Predator = case_when(Specific.cause %in% c("Bobcat", 
                                                    "Coyote", 
                                                    "Lynx", 
                                                    "Avian", 
                                                    "Unknown") ~ Specific.cause,
                              Specific.cause == "Goshawk" ~ "Avian",
                              Specific.cause == "Felid" ~ "Unknown mammal",        # do this for simplicity
                              Specific.cause == "Mammalian" ~ "Unknown mammal"))

# reorder and label factor levels  
fates.pred$Predator <- factor(fates.pred$Predator,
                              levels = c("Lynx",
                                         "Bobcat",
                                         "Unknown felid",
                                         "Coyote",
                                         "Unknown mammal",
                                         "Avian",
                                         "Unknown"))

ggplot(data = fates.pred,
       aes(x = study.year.week,
           fill = Predator)) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "plasma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 7)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:7)

#_______________________________________________________________________________________________
# 4d. Other mortality ----
#_______________________________________________________________________________________________

# reorder and label factor levels  
fates.oth$Cause <- factor(fates.oth$General.cause,
                          levels = c("Unknown",
                                     "Trauma",
                                     "Accident",
                                     "Disease"))

ggplot(data = fates.oth,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "magma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 6)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:6)

#_______________________________________________________________________________________________
# 5. Weekly bar charts - faceted by cluster ----
#_______________________________________________________________________________________________
# 5a. Censoring ----
#_______________________________________________________________________________________________

# factor labels
fates.censor$cluster <- factor(fates.censor$cluster,
                               labels = c("Cluster 1",
                                          "Cluster 2",
                                          "Cluster 3",
                                          "Cluster 4"))

ggplot(data = fates.censor,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # facet
  facet_wrap(~ cluster) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1) +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.24, 2)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:2)

#_______________________________________________________________________________________________
# 5b. Predation mortality ----
#_______________________________________________________________________________________________

# factor labels
fates.pred$cluster <- factor(fates.pred$cluster,
                             labels = c("Cluster 1",
                                        "Cluster 2",
                                        "Cluster 3",
                                        "Cluster 4"))

ggplot(data = fates.pred,
       aes(x = study.year.week,
           fill = Predator)) +
  
  # facet
  facet_wrap(~ cluster) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "plasma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 4)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:4)

#_______________________________________________________________________________________________
# 5c. Other mortality ----
#_______________________________________________________________________________________________

# factor labels
fates.oth$cluster <- factor(fates.oth$cluster,
                            labels = c("Cluster 1",
                                       "Cluster 2",
                                       "Cluster 3",
                                       "Cluster 4"))

ggplot(data = fates.oth,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # facet
  facet_wrap(~ cluster) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "magma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 3)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:3)

#_______________________________________________________________________________________________
# 6. Weekly bar charts - faceted by treatment ----
#_______________________________________________________________________________________________
# 6a. Censoring ----
#_______________________________________________________________________________________________

# factor labels
fates.censor$trt <- factor(fates.censor$trt,
                           labels = c("control",
                                      "retention",
                                      "piling"))

fates.censor$post.trt <- factor(fates.censor$post.trt,
                                labels = c("before",
                                           "after"))

ggplot(data = fates.censor,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # facet
  facet_grid(trt ~ post.trt) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1) +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.24, 2)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:2)

#_______________________________________________________________________________________________
# 6b. Predation mortality ----
#_______________________________________________________________________________________________

# factor labels
fates.pred$trt <- factor(fates.pred$trt,
                         labels = c("control",
                                    "retention",
                                    "piling"))

fates.pred$post.trt <- factor(fates.pred$post.trt,
                              labels = c("before",
                                         "after"))

ggplot(data = fates.pred,
       aes(x = study.year.week,
           fill = Predator)) +
  
  # facet
  facet_grid(trt ~ post.trt) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "plasma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 4)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:4)

#_______________________________________________________________________________________________
# 6c. Other mortality ----
#_______________________________________________________________________________________________

# factor labels
fates.oth$trt <- factor(fates.oth$trt,
                         labels = c("control",
                                    "retention",
                                    "piling"))

fates.oth$post.trt <- factor(fates.oth$post.trt,
                              labels = c("before",
                                         "after"))

ggplot(data = fates.oth,
       aes(x = study.year.week,
           fill = Cause)) +
  
  # facet
  facet_grid(trt ~ post.trt) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks,
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "magma") +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks,
                     labels = month.cutoffs$labels) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 3)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:3)

#_______________________________________________________________________________________________
# 7. All events - faceted by treatment ----
#_______________________________________________________________________________________________
# 7a. Clean data ----
#_______________________________________________________________________________________________

# bind together
fates.all <- rbind(fates.censor[ , 1:16], fates.pred[ , 1:16], fates.oth[ , 1:16])

# add new column separating censors, predation morts, and other morts
fates.all <- fates.all %>%
  
  mutate(Event = case_when(Event.type == "Censor" ~ "Censor",
                           Event.type == "Mortality" & General.cause == "Predation" ~ "Predation",
                           Event.type == "Mortality" & General.cause != "Predation" ~ "Other mortality"))

# reorder factor
fates.all$Event <- factor(fates.all$Event,
                          levels = c("Censor",
                                     "Predation",
                                     "Other mortality"))

#_______________________________________________________________________________________________
# 7b. Weekly ----
#_______________________________________________________________________________________________

ggplot(data = fates.all,
       aes(x = study.year.week,
           fill = Event)) +
  
  # facet
  facet_grid(trt ~ post.trt) +
  
  # white background
  theme_bw() +
  
  # vertical lines at month cutoffs
  geom_vline(xintercept = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
             color = "lightgray") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "mako",
                       begin = 0.3, 
                       end = 1) +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm")
  ) +
  
  # sensible week breaks and labels
  scale_x_continuous(breaks = month.cutoffs$breaks[c(1, 3, 5, 7, 9, 11)],
                     labels = month.cutoffs$labels[c(1, 3, 5, 7, 9, 11)]) +
  
  # coordinates
  coord_cartesian(xlim = c(2.5, 51),
                  ylim = c(0.33, 6)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:6)

#_______________________________________________________________________________________________
# 7c. Monthly ----
#_______________________________________________________________________________________________

fates.all.month <- fates.all %>%
  
  mutate(month = substr(Estimated.event.date, 1, 7)) %>%
  
  # add identifier by cluster groups
  mutate(cluster.group = factor(ifelse(cluster == 4,
                                       "Low",
                                       "High")))

# plot
ggplot(data = fates.all.month,
       aes(x = month,
           fill = cluster.group)) +
  
  # facet
  facet_grid(trt ~ Event) +
  
  # white background
  theme_bw() +
  
  # vertical line for treatment
  geom_vline(xintercept = "2023-10",
             linetype = "dashed") +
  
  # count bars
  geom_bar(stat = "count",
           color = "black") +
  
  # y axis title
  ylab("Number of events") +
  
  # fill label
  labs(fill = "Cluster group") +
  
  # colors
  scale_fill_viridis_d(direction = -1,
                       option = "mako",
                       begin = 0.3, 
                       end = 1) +
  
  # theme arguments
  theme(
    # remove gridlines
    panel.grid = element_blank(),
    
    # adjust axis text
    axis.text.x = element_text(angle = 270,
                               vjust = 0.25),
    
    # remove x axis title
    axis.title.x = element_blank(),
    
    # legend position
    legend.position = "top",
    
    # legend size
    legend.key.size = unit(0.35, "cm"),
    
    # strip text justification
    strip.text = element_text(hjust = 0)
    
  ) +
  
  # sensible week breaks and labels
  scale_x_discrete(breaks = c("2022-10", "2023-02", "2023-06", 
                              "2023-10", "2024-02", "2024-06", 
                              "2024-10"),
                   labels = c("2022-10", "2023-02", "2023-06", 
                              "2023-10", "2024-02", "2024-06", 
                              "2024-10")) +
  
  # coordinates
  coord_cartesian(ylim = c(0.33, 4)) +
  
  # and scale y axis correctly
  scale_y_continuous(breaks = 1:4)
