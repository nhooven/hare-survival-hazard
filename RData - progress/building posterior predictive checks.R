# simulate mortality events 
# we'll use empirical starting points and time-varying covariates from observations

# create blank "fates" datasets
blank.fates <- data.frame()

for (i in 1:length(unique(fates.1$Ear.tag))) {
  
  # subset
  focal.indiv <- unique(fates.1$Ear.tag)[i]
  
  focal.df <- fates.1 %>% filter(Ear.tag == focal.indiv)
  
  # add rows to total to 78 (a year and a half, for buffering purposes)
  # how many rows to add?
  n.rows.add <- 78 - nrow(focal.df)
  
  # weeks recycling
  max.week <- focal.df$week[nrow(focal.df)]
  weeks <- (max.week + 1):((max.week + 1) + n.rows.add - 1)
  weeks.1 <- ifelse(weeks > 52 & weeks <= 104,
                    weeks - 52,
                    ifelse(weeks > 104,
                           weeks - 104,
                           weeks))
  
  # create df to bind (propagate final value)
  add.df <- data.frame(cluster = focal.df$cluster[nrow(focal.df)],
                       Site = focal.df$Site[nrow(focal.df)],
                       Ear.tag = focal.df$Ear.tag[nrow(focal.df)],
                       Collar.type = focal.df$Collar.type[nrow(focal.df)],
                       week = weeks.1,
                       year = NA,
                       status.num = NA,
                       Sex.1 = focal.df$Sex.1[nrow(focal.df)],
                       Mass.1 = focal.df$Mass.1[nrow(focal.df)],
                       HFL.1 = focal.df$HFL.1[nrow(focal.df)],
                       PC1 = focal.df$PC1[nrow(focal.df)],
                       BCI.1 = focal.df$BCI.1[nrow(focal.df)],
                       Treatment.Retention = focal.df$Treatment.Retention[nrow(focal.df)],
                       Treatment.Piling = focal.df$Treatment.Piling[nrow(focal.df)])
  
  # bind together
  focal.df.1 <- rbind(focal.df,
                      add.df)
  
  # remove status.num column
  focal.df.1 <- focal.df.1 %>% dplyr::select(-status.num)
  
  # bind into master df
  blank.fates <- rbind(blank.fates, focal.df.1)
    
}

# simulate survival times using the survival function
# basis functions
knot.list <- quantile(fates.1$week, 
                      probs = seq(from = 0, 
                                  to = 1, 
                                  length.out = 5))

basis.pred <- bs(1:52,       
                 knots = knot.list[-c(1, n.knots)],          
                 degree = 3,                     
                 intercept = FALSE)

# dataset of expected hazard rates for each individual-week combination

# HERE IT MIGHT BE WORTHWHILE TO SUBSAMPLE THE POSTERIOR DRAWS TO SOMETHING EASIER TO
# WORK WITH, SO I CAN RUN THROUGH THE WHOLE PROCESS TO SEE IF IT WORKS
model.draws.sample <- model.draws %>%
  
  slice(sample(1:4000,
               size = 1000))

all.draws.indivs <- data.frame()

# time it
start.time <- Sys.time()

# loop
for (i in 1:nrow(model.draws.sample)) {
  
  focal.draw <- model.draws.sample[i, ]
  
  # h0(t) (baseline hazard - spline)
  # spline weights
  w0 <- t(as.matrix(as.numeric(focal.draw[ , 7:12])))
  
  # multiply with 'sweep'
  w0.by.b <- sweep(basis.pred, 2, w0, `*`)
  
  # sum to create "normal" scale spline prediction
  w0.by.b.sum <- apply(w0.by.b, 1, sum)
  
  # calculate for every observation
  blank.fates.1 <- 
    
    mutate(blank.fates, 
           
           lambda = 
             
             # calculate total expected hazard
             as.numeric(exp((
               
                             # base intercept
                             focal.draw$a0 + 
                               
                             # standard deviation
                             focal.draw$sigma * 
                             
                             # non-centered random intercept term
                             focal.draw[2:5][blank.fates$cluster]) + 
                             
                             # extract correct weekly hazard ("normal" scale)
                             w0.by.b.sum[blank.fates$week])) * 
               
               # covariate effects 
               exp(blank.fates$Sex.1 * focal.draw$b_sex +
                   blank.fates$Mass.1 * focal.draw$b_mas +
                   blank.fates$HFL.1 * focal.draw$b_hfl +
                   blank.fates$BCI.1 * focal.draw$b_bci +
                   blank.fates$Treatment.Retention * focal.draw$b_ret +
                   blank.fates$Treatment.Piling * focal.draw$b_pil),
           
           # add draw number
           draw = i)
  
  # calculate cumulative hazard, survival, and fit function
  # nested loop (by individual)
  all.indivs <- data.frame()
  
  for (j in 1:length(unique(blank.fates.1$Ear.tag))) {
    
    focal.indiv <- unique(blank.fates.1$Ear.tag)[j]
    
    focal.indiv.df <- blank.fates.1 %>% filter(Ear.tag == focal.indiv)
    
    # add cumulative sum
    focal.indiv.df <- focal.indiv.df %>% 
      
      mutate(cumul.haz = cumsum(lambda)) %>%
      
      # convert to cumulative survival and add sequential week
      mutate(cumul.surv = exp(-cumul.haz),
             seq.week = as.numeric(rownames(focal.indiv.df)))
    
    # fit smooth function to cumulative survival curve
    surv.spline <- lm(data = focal.indiv.df,
                      formula = cumul.surv ~ bs(x = seq.week,
                                                df = 10))      # this provides decent flexibility
    
    # sample from a uniform distribution, snap to nearest predicted y/extract the x
    sampled.y <- runif(n = 1, 
                       min(predict(surv.spline)), 
                       max(predict(surv.spline)))
    
    snap <- as.integer(which.min(abs(predict(surv.spline) - sampled.y)))
    
    # bind into data frame with characteristics of first observation
    focal.indiv.df.1 <- focal.indiv.df %>%
      
      # slice the first row
      slice(1) %>%
      
      # keep only columns we need
      dplyr::select(-c(lambda,
                       cumul.haz,
                       cumul.surv,
                       seq.week)) %>%
      
      # add in lifetime (in weeks)
      mutate(lifetime = snap)
    
    # bind into sub-master df
    all.indivs <- rbind(all.indivs, focal.indiv.df.1)
      
  }
  
  # bind into master df
  all.draws.indivs <- rbind(all.draws.indivs, all.indivs)
  
}

# end time
Sys.time() - start.time

# plot
ggplot(data = all.draws.indivs) +
  
  theme_classic() +
  
  geom_line(aes(x = lifetime,
                   group = draw),
               stat = "density",
               color = "#33CCCC",
               alpha = 0.15) +
  
  # empirical distribution
  geom_density(data = lifetimes.empirical,
               aes(x = lifetime),
               color = "black",
               fill = NA,
               linewidth = 1.15)

# this is clearly leading to FAR too many individuals surviving a long time,
# especially > 30 weeks
# we must incorporate the censoring process here!
# either in the model (eventually), or in the data thinning here
