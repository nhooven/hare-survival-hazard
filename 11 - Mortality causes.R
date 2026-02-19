# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 11 - Mortality causes & data summary plots
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 18 Feb 2026
# Date completed: 19 Feb 2026
# Date last modified: 19 Feb 2026
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data

#_______________________________________________________________________________________________
# 2. Read in fates data ----
#_______________________________________________________________________________________________

fates <- read.csv("Raw data/fates_02_18_2026.csv")

#_______________________________________________________________________________________________
# 3. Clean data ----
#_______________________________________________________________________________________________

fates.1 <- fates %>%
  
  # keep relevant columns
  dplyr::select(Site,
                Sex,
                MRID,
                Collar.type,
                Status,
                Capture.date,
                Estimated.event.date,
                Event.type,
                General.cause,
                Specific.cause,
                Include) %>%
  
  # filter mortalities
  filter(Event.type == "Mortality") %>%
  
  # parse dates
  mutate(cap.date = mdy(Capture.date),
         mort.date = mdy(Estimated.event.date)) %>%
  
  # remove original date columns
  dplyr::select(-Capture.date, -Estimated.event.date) %>%
  
  # cluster and treatment variables
  mutate(
    
    cluster = substr(Site, 1, 1),
    
    trt = case_when(Site %in% c("1C", "2C", "3C", "4C") ~ "untreated",
                    Site %in% c("1A", "2B", "3B", "4A") ~ "retention",
                    Site %in% c("1B", "2A", "3A", "4B") ~ "piling")
    
    ) %>%
  
  # post-treatment
  mutate(
    
    post1 = case_when(
      
      cluster == 1 & mort.date > ymd("2023-10-12") & mort.date < ymd("2024-10-12") ~ 1,
      cluster == 2 & mort.date > ymd("2023-10-08") & mort.date < ymd("2024-10-08") ~ 1,
      cluster == 3 & mort.date > ymd("2023-10-05") & mort.date < ymd("2024-10-05") ~ 1,
      cluster == 4 & mort.date > ymd("2023-10-10") & mort.date < ymd("2024-10-10") ~ 1),
    
    post2 = case_when(
      
      cluster == 1 & mort.date >= ymd("2024-10-12") ~ 1,
      cluster == 2 & mort.date >= ymd("2024-10-08") ~ 1,
      cluster == 3 & mort.date >= ymd("2024-10-05") ~ 1,
      cluster == 4 & mort.date >= ymd("2024-10-10") ~ 1)
    
  ) %>%
  
  # switch NAs to 0
  replace_na(list(post1 = 0,
                  post2 = 0))

#_______________________________________________________________________________________________
# 4. Tabulate ----
#_______________________________________________________________________________________________
# 4a. Overall individuals collared ----
#_______________________________________________________________________________________________

all.indivs <- fates %>% 
  
  group_by(Animal.ID) %>%
  
  summarize(n())

# sex
fates %>%
  
  group_by(Animal.ID) %>%
  
  slice(1) %>%
  
  ungroup %>% 
  
  group_by(Sex) %>%
  
  summarize(n())

#_______________________________________________________________________________________________
# 4b. Overall mortality causes ----
#_______________________________________________________________________________________________

overall.mort <- fates.1 %>%
  
  group_by(General.cause) %>%
  
  summarize(n())

overall.mort

sum(overall.mort$`n()`)

#_______________________________________________________________________________________________
# 4c. Mortality causes (for model) ----
#_______________________________________________________________________________________________

# overall
kept.mort <- fates.1 %>%

  filter(Include == "Y") %>%
  
  group_by(General.cause) %>%
  
  summarize(n())

kept.mort

sum(kept.mort$`n()`)

#_______________________________________________________________________________________________
# 4d. Uninformative censoring ----
#_______________________________________________________________________________________________

uninf.cens <- fates %>%
  
  filter(Event.type == "Censor" &
         General.cause == "Unknown" &
         Include == "Y")

#_______________________________________________________________________________________________
# 5. Specific causes, including predator identity ----
#_______________________________________________________________________________________________

fates.2 <- fates.1 %>%
  
  filter(Include == "Y") %>%
  
  # replace blank specific cause with "unknown"
  mutate(Specific.cause = ifelse(Specific.cause == "", "Unknown", Specific.cause)) %>%
  
  # life zone
  mutate(lz = factor(ifelse(cluster == 4, "XMC", "SFL"),
                     labels = c("a", "b"))) %>%
  
  # levels of classification
  mutate(
    
    # CLASS
    cause.class = case_when(
      
      # non-predation
      General.cause == "Unknown" ~ "Unknown",
      Specific.cause == "Covered in snow" ~ "Accident",
      Specific.cause == "Infection" ~ "Disease",
      Specific.cause == "Not predation" ~ "Unknown",
      General.cause == "Trauma" ~ "Trauma",
      
      # predation
      General.cause == "Predation" & Specific.cause == "Unknown" ~ "Unknown predation",
      Specific.cause %in% c("Bobcat", "Coyote", "Felid", "Lynx", "Mammalian") ~ "Mammalian",
      Specific.cause %in% c("Avian", "Goshawk") ~ "Avian"
      
    ),
    
    # FAMILY
    cause.family = case_when(
      
      # non-predation
      General.cause == "Unknown" ~ "Unknown",
      Specific.cause == "Covered in snow" ~ "Accident",
      Specific.cause == "Infection" ~ "Disease",
      Specific.cause == "Not predation" ~ "Unknown",
      General.cause == "Trauma" ~ "Trauma",
      
      # predation
      General.cause == "Predation" & Specific.cause == "Unknown" ~ "Unknown predation",
      Specific.cause %in% c("Avian", "Goshawk") ~ "Avian",
      Specific.cause %in% c("Bobcat", "Felid", "Lynx") ~ "Felid",
      Specific.cause == "Coyote" ~ "Canid",
      Specific.cause == "Mammalian" ~ "Unknown mammal"
    
    ),
    
    # SPECIES
    cause.species = case_when(
      
      
      # non-predation
      General.cause == "Unknown" ~ "Unknown",
      Specific.cause == "Covered in snow" ~ "Accident",
      Specific.cause == "Infection" ~ "Disease",
      Specific.cause == "Not predation" ~ "Unknown",
      General.cause == "Trauma" ~ "Trauma",
      
      # predation
      General.cause == "Predation" & Specific.cause == "Unknown" ~ "Unknown predation",
      Specific.cause %in% c("Avian", "Goshawk") ~ "Avian",
      Specific.cause == "Felid" ~ "Unknown felid",
      Specific.cause == "Bobcat" ~ "Bobcat",
      Specific.cause == "Lynx" ~ "Lynx",
      Specific.cause == "Coyote" ~ "Coyote",
      Specific.cause == "Mammalian" ~ "Unknown mammal"
    
   )
  
 ) %>%
  
  # reorder general cause factor
  mutate(General.cause = factor(General.cause,
                                levels = c("Predation",
                                           "Unknown",
                                           "Trauma",
                                           "Accident",
                                           "Disease"))) %>%
  
  # add month
  mutate(month = substr(mort.date, 6, 7)) %>%
  mutate(month = factor(month,
                        levels = c("10", "11", "12", "01", "02", "03", "04", "05", "06", "07", "08", "09"),
                        labels = c("Oct", 
                                   "Nov", 
                                   "Dec",
                                   "Jan", 
                                   "Feb", 
                                   "Mar", 
                                   "Apr", 
                                   "May", 
                                   "Jun", 
                                   "Jul", 
                                   "Aug", 
                                   "Sep")))

# summarize
fates.2 %>%
  
  filter(General.cause == "Predation") %>%
  
  group_by(cause.class, cause.family, cause.species) %>%
  
  summarize(n())

# n predator species
1 + 15 + 9 + 39

#_______________________________________________________________________________________________
# 6. Count plots ----
#_______________________________________________________________________________________________
# 6a. Predation vs. all others, split by life zone ----
#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ lz, nrow = 2) +
  
  geom_bar(data = fates.2,
           aes(x = month,
               fill = General.cause)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  
  scale_fill_viridis_d() +
  
  # x axis labels
  scale_x_discrete(labels = c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S"))
  
#_______________________________________________________________________________________________
# 6b. Predation vs. all others, treatment ----

# factor variables for facetting
fates.3 <- fates.2 %>%
  
  mutate(
    
    post = factor(
      
      case_when(
        
        post1 == 1 ~ "POST1",
        post2 == 1 ~ "POST2",
        post1 == 0 & post2 == 0 ~ "PRE"
        
        ),
      
      levels = c("PRE", "POST1", "POST2")
      
      ),
    
    trt = factor(trt, levels = c("untreated", "retention", "piling"))
  
  )

#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_grid(trt ~ post) +
  
  geom_bar(data = fates.3,
           aes(x = month,
               fill = General.cause)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  
  scale_fill_viridis_d() +
  
  xlab("Month of year") +
  
  # x axis labels
  scale_x_discrete(labels = c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S"))

#_______________________________________________________________________________________________
# 6c. Predator identity, treatment ----

# new variable
fates.4 <- fates.3 %>%
  
  mutate(
    
    cause.species.1 = factor(
      
      ifelse(
        
        General.cause %in% c("Unknown", "Trauma", "Accident", "Disease"),
        "Other",
        cause.species
        
        ),
      
      levels = c("Lynx", 
                 "Coyote", 
                 "Bobcat", 
                 "Unknown felid", 
                 "Avian", 
                 "Unknown mammal", 
                 "Unknown predation",
                 "Other")
      
      )
  
  )

#_______________________________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_grid(trt ~ post) +
  
  geom_bar(data = fates.4,
           aes(x = month,
               fill = cause.species.1)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 7),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  
  scale_fill_manual(values = c("#000004", "#320a5e", "#781c6d", "#bc3754", "#ed6925", "#fbb61a", "#fcffa4", "lightgray")) +
  
  xlab("Month of year") +
  ylab("Mortalities") +
  
  # x axis labels
  scale_x_discrete(labels = c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S")) +
  
  scale_y_continuous(breaks = c(1, 3, 5, 7)) +
  
  coord_cartesian(ylim = c(0.4, 7))

# 823 x 423

#_______________________________________________________________________________________________
# 7. Monitoring plots ----

# this will depict our sample size across the entire study

#_______________________________________________________________________________________________
# 7a. Set up data frame ----
#_______________________________________________________________________________________________

monitoring <- fates %>%
  
  # keep relevant columns
  dplyr::select(Site,
                Sex,
                MRID,
                Collar.type,
                Status,
                Capture.date,
                Estimated.event.date,
                Event.type,
                General.cause,
                Specific.cause,
                Include) %>%
  
  # filter only those included
  filter(Include == "Y") %>%
  
  # parse dates
  mutate(cap.date = mdy(Capture.date),
         mort.date = mdy(Estimated.event.date)) %>%
  
  # remove original date columns
  dplyr::select(-Capture.date, -Estimated.event.date) %>%
  
  # deployment, cluster, and treatment variables
  mutate(
    
    deployment = paste0(MRID, "_", cap.date),
    
    cluster = substr(Site, 1, 1),
    
    trt = factor(case_when(Site %in% c("1C", "2C", "3C", "4C") ~ "untreated",
                           Site %in% c("1A", "2B", "3B", "4A") ~ "retention",
                           Site %in% c("1B", "2A", "3A", "4B") ~ "piling"),
                 
                 levels = c("untreated", "retention", "piling"))
    
  )

#_______________________________________________________________________________________________
# 7b. Full plot ----

# we'll stack these instead of facetting
library(cowplot)

#_______________________________________________________________________________________________

# untreated
monitor.1 <- ggplot(data = monitoring %>% filter(trt == "untreated")) +
  
  theme_bw() +
  
  # add year boundaries
  geom_vline(xintercept = c(ymd("2023-01-01"),
                            ymd("2024-01-01"),
                            ymd("2025-01-01")),
             linetype = "dashed",
             color = "gray65") +
  
  # add treatment
  geom_rect(xmin = ymd("2023-10-05"),
            xmax = ymd("2023-10-12"),
            ymin = -4,
            ymax = 384,
            fill = "gray75") +
  
  # add segments
  geom_segment(aes(x = cap.date,
                   xend = mort.date,
                   y = reorder(deployment, 
                               cap.date,
                               decreasing = T),
                   color = Event.type)) +
  
  scale_color_manual(values = c("gray50", "#FF3300")) +
  
  # add end points
  geom_point(aes(x = mort.date,
                 y = reorder(deployment, 
                             cap.date,
                             decreasing = T),
                 shape = Event.type),
             size = 0.8,
             color = "#FF3300",
             stroke = 0.8) +
  
  scale_shape_manual(values = c(NA, 4)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.96, vjust = -10),
        plot.margin = margin(t = -15, b = 0)) +
  
  # axis padding
  coord_cartesian(xlim = c(ymd("2022-11-15"), ymd("2025-09-15")),
                  ylim = c(-3, 118)) +
  
  ggtitle("a")

# retention
monitor.2 <- ggplot(data = monitoring %>% filter(trt == "retention")) +
  
  theme_bw() +
  
  # add year boundaries
  geom_vline(xintercept = c(ymd("2023-01-01"),
                            ymd("2024-01-01"),
                            ymd("2025-01-01")),
             linetype = "dashed",
             color = "gray65") +
  
  # add treatment
  geom_rect(xmin = ymd("2023-10-05"),
            xmax = ymd("2023-10-12"),
            ymin = -4,
            ymax = 384,
            fill = "gray75") +
  
  # add segments
  geom_segment(aes(x = cap.date,
                   xend = mort.date,
                   y = reorder(deployment, 
                               cap.date,
                               decreasing = T),
                   color = Event.type)) +
  
  scale_color_manual(values = c("gray50", "#FF3300")) +
  
  # add end points
  geom_point(aes(x = mort.date,
                 y = reorder(deployment, 
                             cap.date,
                             decreasing = T),
                 shape = Event.type),
             size = 0.8,
             color = "#FF3300",
             stroke = 0.8) +
  
  scale_shape_manual(values = c(NA, 4)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.96, vjust = -10),
        plot.margin = margin(t = -15, b = 0)) +
  
  # axis padding
  coord_cartesian(xlim = c(ymd("2022-11-15"), ymd("2025-09-15")),
                  ylim = c(-3, 145)) +
  
  ggtitle("b")

# piling
monitor.3 <- ggplot(data = monitoring %>% filter(trt == "piling")) +
  
  theme_bw() +
  
  # add year boundaries
  geom_vline(xintercept = c(ymd("2023-01-01"),
                            ymd("2024-01-01"),
                            ymd("2025-01-01")),
             linetype = "dashed",
             color = "gray65") +
  
  # add treatment
  geom_rect(xmin = ymd("2023-10-05"),
            xmax = ymd("2023-10-12"),
            ymin = -4,
            ymax = 384,
            fill = "gray75") +
  
  # add segments
  geom_segment(aes(x = cap.date,
                   xend = mort.date,
                   y = reorder(deployment, 
                               cap.date,
                               decreasing = T),
                   color = Event.type)) +
  
  scale_color_manual(values = c("gray50", "#FF3300")) +
  
  # add end points
  geom_point(aes(x = mort.date,
                 y = reorder(deployment, 
                             cap.date,
                             decreasing = T),
                 shape = Event.type),
             size = 0.8,
             color = "#FF3300",
             stroke = 0.8) +
  
  scale_shape_manual(values = c(NA, 4)) +
  
  # theme
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.96, vjust = -10),
        plot.margin = margin(t = -15, b = 0)) +
  
  # x axis breaks - middle of years
  scale_x_date(breaks = c(ymd("2022-11-01"), ymd("2023-06-30"), ymd("2024-06-30"), ymd("2025-07-10")),
               labels = c("2022", "2023", "2024", "2025")) +
  
  # axis padding
  coord_cartesian(xlim = c(ymd("2022-11-15"), ymd("2025-09-15")),
                  ylim = c(-3, 132)) +
  
  ggtitle("c")

# plot together
plot_grid(monitor.1, monitor.2, monitor.3, nrow = 3, rel_heights = c(1, 1, 1.09))

# 465 x 735

#_______________________________________________________________________________________________
# 8. Cause-specific mortality table ----
#_______________________________________________________________________________________________

csm.table <- fates.3 %>%
  
  group_by(General.cause,
           cause.class,
           cause.family,
           cause.species,
           trt,
           post) %>%
  
  summarize(n()) %>%
  
  ungroup() %>%
  
  pivot_wider(names_from = c(trt, post),
              values_from = `n()`)

csm.table

write.table(csm.table, "clipboard", sep = "\t")
