# Project: WSU Snowshoe Hare and PCT Project
# Subproject: Survival and hazard modeling
# Script: 09 - Visualization - Baseline hazard
# Author: Nathan D. Hooven, Graduate Research Assistant
# Email: nathan.hooven@wsu.edu / nathan.d.hooven@gmail.com
# Date began: 17 Nov 2025 
# Date completed: 02 Mar 2026 
# Date last modified: 02 Mar 2026 
# R version: 4.4.3

#_______________________________________________________________________________________________
# 1. Load required packages ----
#_______________________________________________________________________________________________

library(tidyverse)       # manipulate and clean data
library(mgcv)

#_______________________________________________________________________________________________
# 2. Read in and clean data ----
#_______________________________________________________________________________________________

# model samples
model.fit.1 <- readRDS("Model outputs/model_1.rds")
model.fit.2 <- readRDS("Model outputs/model_2.rds")
model.fit.3 <- readRDS("Model outputs/model_3.rds")

# dataset
fates <- read.csv("Cleaned data/fates_forModel.csv")

# convert to data.frame
model.fit.1 <- as.data.frame(do.call(rbind, model.fit.1))
model.fit.2 <- as.data.frame(do.call(rbind, model.fit.2))
model.fit.3 <- as.data.frame(do.call(rbind, model.fit.3))

#_______________________________________________________________________________________________
# 3. Define spline ----
#_______________________________________________________________________________________________

# define knots (quantile)
n.knots <- 9 + 1

knot.list <- quantile(fates$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = n.knots))

# weeks to predict on
weeks.pred <- seq(1, 52, length.out = 100)

# first, we need to create the "normal" scale predictions
# let's create a prediction df first
basis.pred <- cSplineDes(weeks.pred,
                         knots = knot.list,
                         ord = 4)

#_______________________________________________________________________________________________
# 4. Calculate each hazard prediction ----
#_______________________________________________________________________________________________

# function (accepts model fit)
haz_spline <- function (
    
  x,
  low = 0.025,
  upp = 0.975
  
  ) {
  
  # keep only columns we need
  x.1 <- x %>%
    
    dplyr::select(c("a0sc[1, 1]":"a0sc[2, 4]",
                    "wsc[1, 1, 1]":"wsc[2, 4, 9]"))
  
  # calculate predictions for each sex-forest type group
  # internal function to calculate predictions
  spline_pred <- function (
    
    y, 
    sex,         # 1 == F, 2 == M
    cluster
    
    ) {
    
    # define a blank matrix [n.draws, n.prediction points]
    pred.matrix <- matrix(data = NA,
                          nrow = nrow(x.1),
                          ncol = nrow(basis.pred))
    
    # convert coefficients to column vector for multiplication
    w <- matrix(data = c(y[[paste0("wsc[", sex, ", ", cluster, ", 1]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 2]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 3]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 4]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 5]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 6]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 7]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 8]")]],
                         y[[paste0("wsc[", sex, ", ", cluster, ", 9]")]]),
           ncol = 9)

    
    # multiply with 'sweep'
    w.by.b <- sweep(basis.pred, 2, w, `*`)
    
    # sum to create "normal" scale spline prediction
    w.by.b.sum <- apply(w.by.b, 1, sum) 
    
    # add to the intercept (a0) and exponentiate to transform to hazard rate scale
    spline.pred <- exp(as.numeric(y[paste0("a0sc[", sex, ", ", cluster, "]")]) + w.by.b.sum) 
    
    # and add into the array
    return(spline.pred)
    
  }
  
  # calculate splines
  spl.1.1 <- apply(x.1, 1, spline_pred, 1, 1)
  spl.1.2 <- apply(x.1, 1, spline_pred, 1, 2)
  spl.1.3 <- apply(x.1, 1, spline_pred, 1, 3)
  spl.1.4 <- apply(x.1, 1, spline_pred, 1, 4)
  spl.2.1 <- apply(x.1, 1, spline_pred, 2, 1)
  spl.2.2 <- apply(x.1, 1, spline_pred, 2, 2)
  spl.2.3 <- apply(x.1, 1, spline_pred, 2, 3)
  spl.2.4 <- apply(x.1, 1, spline_pred, 2, 4)
  
  # summarize function
  summ_spline <- function (z, sex, cluster) {
    
    return(
      
      data.frame(
        
        med = apply(z, 1, median),
        low = apply(z, 1, quantile, prob = low),
        upp = apply(z, 1, quantile, prob = upp),
        sex = sex,
        cluster = cluster
      
        )
    
    )
    
  }
  
  # bind together summaries
  all.summaries <- rbind(
    
    summ_spline(spl.1.1, 1, 1),
    summ_spline(spl.1.2, 1, 2),
    summ_spline(spl.1.3, 1, 3),
    summ_spline(spl.1.4, 1, 4),
    summ_spline(spl.2.1, 2, 1),
    summ_spline(spl.2.2, 2, 2),
    summ_spline(spl.2.3, 2, 3),
    summ_spline(spl.2.4, 2, 4)
    
  )
  
  # calculate hazard ratios FIRST on the full posterior predictions
  hr.1 <- spl.2.1 / spl.1.1
  hr.2 <- spl.2.2 / spl.1.2
  hr.3 <- spl.2.3 / spl.1.3
  hr.4 <- spl.2.4 / spl.1.4
  
  # summarize HRs
  hr.summaries <- rbind(
    
    summ_spline(hr.1, "HR", 1),
    summ_spline(hr.2, "HR", 2),
    summ_spline(hr.3, "HR", 3),
    summ_spline(hr.4, "HR", 4)
    
  )
  
  # make a list
  spline.list <- list(all.summaries,
                      hr.summaries)
  
  # return
  return(spline.list)
  
}

# apply function
spline.summary.mod1 <- haz_spline(model.fit.1)
spline.summary.mod2 <- haz_spline(model.fit.2)
spline.summary.mod3 <- haz_spline(model.fit.3)

#_______________________________________________________________________________________________
# 4. Process predictions for plotting ----
#_______________________________________________________________________________________________
# 4a. Sex x cluster ----
#_______________________________________________________________________________________________

# function
spline_forPlot <- function (x) {
  
  x.1 <- x %>%
    
    mutate(
      
      # add week
      t = rep(weeks.pred, times = 8),
      
      # add life zone
      lz = ifelse(cluster == 4, "XMC", "SFL"),
      
      # change labels
      sex = factor(sex, labels = c("sex == F", "sex == M")),
      cluster = factor(cluster, labels = c("Cluster 1",
                                           "Cluster 2",
                                           "Cluster 3",
                                           "Cluster 4"))
      
    )
  
  return(x.1)
  
}

# apply function
spline.forPlot.mod1 <- spline_forPlot(spline.summary.mod1[[1]])
spline.forPlot.mod2 <- spline_forPlot(spline.summary.mod2[[1]])
spline.forPlot.mod3 <- spline_forPlot(spline.summary.mod3[[1]])

#_______________________________________________________________________________________________
# 4b. M:F hazard ratio ----
#_______________________________________________________________________________________________

# function
spline_hr_forPlot <- function (x) {
  
  x.1 <- x %>%
    
    mutate(
      
      # add week
      t = rep(weeks.pred, times = 4),
      
      # add life zone
      lz = ifelse(cluster == 4, "XMC", "SFL"),
      
      # change labels
      cluster = factor(cluster, labels = c("Cluster 1",
                                           "Cluster 2",
                                           "Cluster 3",
                                           "Cluster 4"))
      
    )
  
  return(x.1)
  
}

# apply function
spline.hr.forPlot.mod1 <- spline_hr_forPlot(spline.summary.mod1[[2]])
spline.hr.forPlot.mod2 <- spline_hr_forPlot(spline.summary.mod2[[2]])
spline.hr.forPlot.mod3 <- spline_hr_forPlot(spline.summary.mod3[[2]])

#_______________________________________________________________________________________________
# 5. Plot predictions ----
#_______________________________________________________________________________________________

# define month cutoffs
month.cutoffs.all <- data.frame(breaks = c(1, 5.57, 9.86, 
                                       14.29, 18.71, 22.71, 
                                       27.14, 31.43, 35.86, 
                                       40.14, 44.57, 48.86),
                            labels = c("O", "N", "D", 
                                       "J", "F", "M", 
                                       "A", "M", "J", 
                                       "J", "A", "S"))

month.cutoffs.some <- data.frame(breaks = c(5.57, 14.29, 22.71, 31.43, 40.14, 48.86),
                                labels = c("Nov", "Jan", "Mar", "May", "Jul", "Sep"))

#_______________________________________________________________________________________________
# 5a. Expected hazard across the annual cycle, one model at a time ----
#_______________________________________________________________________________________________

# function
haz_plot <- function (x) {
  
  ggplot(data = x) +
    
    theme_bw() +
    
    facet_grid(sex ~ cluster) +
    
    # credible interval
    geom_ribbon(aes(x = t,
                    ymin = low,
                    ymax = upp,
                    fill = lz),
                alpha = 0.25) +
    
    # median
    geom_line(aes(x = t,
                  y = med,
                  color = lz),
              linewidth = 0.9) +
    
    # theme arguments
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(vjust = 0.25,
                                     hjust = 0,
                                     angle = 90,
                                     size = 7),
          legend.position = "none",
          legend.background = element_rect(color = "black"),
          strip.text = element_text(hjust = 0),
          strip.text.x = element_text(size = 7),
          strip.background = element_rect(fill = "white")) +
    
    # axis labels
    scale_x_continuous(breaks = month.cutoffs.some$breaks,
                       labels = month.cutoffs.some$labels) +
    
    scale_y_continuous(breaks = c(0.04, 0.12, 0.20, 0.28)) +
    
    # axis labels
    ylab("Baseline hazard") +
    
    # coordinates
    coord_cartesian(xlim = c(3.5, 49.7),
                    ylim = c(0.009, 0.32)) +
    
    # colors
    scale_color_manual(values = c("#003300", "#669900")) +
    scale_fill_manual(values = c("#003300", "#669900"))
    
}

# apply function
haz_plot(spline.forPlot.mod2)

# 600 x 320

#_______________________________________________________________________________________________
# 5b. Expected hazard across the annual cycle, all models ----
#_______________________________________________________________________________________________

spline.forPlot.all <- rbind(spline.forPlot.mod1 %>% mutate(scenario = "1"),
                            spline.forPlot.mod2 %>% mutate(scenario = "2"),
                            spline.forPlot.mod3 %>% mutate(scenario = "3"))


# function
haz_plot_all <- function (x) {
  
  ggplot(data = x) +
    
    theme_bw() +
    
    facet_grid(sex ~ cluster) +
    
    # credible interval
    geom_ribbon(aes(x = t,
                    ymin = low,
                    ymax = upp,
                    fill = scenario),
                alpha = 0.15) +
    
    # median
    geom_line(aes(x = t,
                  y = med,
                  color = scenario),
              linewidth = 0.9) +
    
    # theme arguments
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(vjust = 0.25,
                                     hjust = 0,
                                     angle = 90,
                                     size = 7),
          legend.position = "top",
          legend.background = element_rect(color = "black"),
          legend.title = element_text(size = 8),
          strip.text = element_text(hjust = 0),
          strip.text.x = element_text(size = 7),
          strip.background = element_rect(fill = "white")) +
    
    # axis labels
    scale_x_continuous(breaks = month.cutoffs.some$breaks,
                       labels = month.cutoffs.some$labels) +
    
    scale_y_continuous(breaks = c(0.04, 0.12, 0.20, 0.28)) +
    
    # axis labels
    ylab("Baseline hazard") +
    
    # coordinates
    coord_cartesian(xlim = c(3.5, 49.7),
                    ylim = c(0.009, 0.32)) +
    
    scale_color_viridis_d() +
    scale_fill_viridis_d()
  
}

# 600 x 360 

#_______________________________________________________________________________________________
# 5b. log(M:F hazard ratio) across the annual cycle ----

# the log scale is really necessary for visualization here
# for scenario 2, we can tighten the y axis

#_______________________________________________________________________________________________

# function
hr_plot <- function (x) {
  
  ggplot(data = x) +
    
    theme_bw() +
    
    facet_wrap(~ cluster, nrow = 1) +
    
    # zero line
    geom_hline(yintercept = 0,
               color = "gray50",
               linetype = "dashed") +
    
    # credible interval
    geom_ribbon(aes(x = t,
                    ymin = log(low),
                    ymax = log(upp),
                    fill = lz),
                alpha = 0.25) +
    
    # median
    geom_line(aes(x = t,
                  y = log(med),
                  color = lz),
              linewidth = 0.9) +
    
    # theme arguments
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(vjust = 0.25,
                                     hjust = 0,
                                     angle = 90,
                                     size = 7),
          legend.position = "none",
          legend.background = element_rect(color = "black"),
          strip.text = element_text(hjust = 0),
          strip.text.x = element_text(size = 7),
          strip.background = element_rect(fill = "white")) +
    
    # axis labels
    scale_x_continuous(breaks = month.cutoffs.some$breaks,
                       labels = month.cutoffs.some$labels) +
    
    scale_y_continuous(breaks = seq(-3, 3, 1)) +
    
    coord_cartesian(xlim = c(3.5, 49.7),
                    ylim = c(-3.5, 2.65)) +
    
    # axis labels
    ylab("ln(hazard ratio)") +
    
    # colors
    scale_color_manual(values = c("#003300", "#669900")) +
    scale_fill_manual(values = c("#003300", "#669900"))
  
}

hr_plot(spline.hr.forPlot.mod2)

# 600 x 200
