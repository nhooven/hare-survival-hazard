# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 11 - Mortality causes
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 18 Feb 2026
# Date completed: 
# Date last modified: 18 Feb 2026
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
# 5. Predator identity ----
#_______________________________________________________________________________________________

fates.pred <- fates.1 %>%
  
  filter(Include == "Y" &
         General.cause == "Predation") %>%
  
  # replace blank specific cause with "unknown"
  mutate(Specific.cause = ifelse(Specific.cause == "", "Unknown", Specific.cause)) %>%
  
  # life zone
  mutate(lz = ifelse(cluster == 4, "XMC", "SFL")) %>%
  
  # levels of classification
  mutate(
    
    # CLASS
    pred.class = case_when(
      
      Specific.cause == "Unknown" ~ "Unknown",
      Specific.cause %in% c("Bobcat", "Coyote", "Felid", "Lynx", "Mammalian") ~ "Mammalian",
      Specific.cause %in% c("Avian", "Goshawk") ~ "Avian"
      
    ),
    
    # FAMILY
    pred.family = case_when(
      
    Specific.cause == "Unknown" ~ "Unknown",
    Specific.cause %in% c("Avian", "Goshawk") ~ "Avian",
    Specific.cause %in% c("Bobcat", "Felid", "Lynx") ~ "Felid",
    Specific.cause == "Coyote" ~ "Canid",
    Specific.cause == "Mammalian" ~ "Unknown mammal"
    
    ),
    
    # SPECIES
    pred.species = case_when(
      
    Specific.cause == "Unknown" ~ "Unknown",
    Specific.cause %in% c("Avian", "Goshawk") ~ "Avian",
    Specific.cause == "Felid" ~ "Unknown felid",
    Specific.cause == "Bobcat" ~ "Bobcat",
    Specific.cause == "Lynx" ~ "Lynx",
    Specific.cause == "Coyote" ~ "Coyote",
    Specific.cause == "Mammalian" ~ "Unknown mammal"
    
  )
  
  )

# summarize
# 
fates.pred %>%
  
  group_by(pred.class, pred.family, pred.species) %>%
  
  summarize(n())

# n predator species
1 + 15 + 9 + 39

#_______________________________________________________________________________________________
# 6. Count plots ----
#_______________________________________________________________________________________________
# 6a. Predation vs. all others, split by life zone ----
#_______________________________________________________________________________________________

# new df for plotting
fates.2 <- fates.1 %>%
  
  filter(Include == "Y") %>%
  
  mutate(month = substr(mort.date, 6, 7),
         lz = ifelse(cluster == 4, "XMC", "SFL")) %>%
  
  # reorder general cause factor
  mutate(General.cause = factor(General.cause,
                                levels = c("Predation",
                                           "Unknown",
                                           "Trauma",
                                           "Accident",
                                           "Disease")))

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ lz, nrow = 2) +
  
  geom_bar(data = fates.2,
           aes(x = month,
               fill = General.cause)) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black")) +
  
  scale_fill_viridis_d()
