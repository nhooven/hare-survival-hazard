# fit the same thing in mgcv

ml.model <- gam(data = fates,
                formula = mort.pred ~ s(week, 
                                        k = 7,
                                        bs = "cc",
                                        fx = FALSE,
                                        m = 1),
                knots = as.list(unname(knot.list)),
                family = poisson())

summary(ml.model)

# spline parameters
ml.model$coefficients
exp(ml.model$coefficients)

# smoothing parameter
ml.model$sp

# predict
ml.pred <- data.frame(week = seq(1, 52, length.out = 1000))

ml.pred$fit <- predict(ml.model, 
                       newdata = ml.pred,
                       se.fit = TRUE)$fit

ml.pred$se.fit <- predict(ml.model, 
                          newdata = ml.pred,
                          se.fit = TRUE)$se.fit

ggplot(data = ml.pred,
       aes(x = week,
           y = exp(fit))) +
  
  geom_line(color = "black") +
  
  geom_ribbon(aes(x = week,
                  y = exp(fit),
                  ymin = exp(fit - se.fit * 1.96),
                  ymax = exp(fit + se.fit * 1.96)),
              alpha = 0.25)


smooth.test.1 <- smooth.construct(s(week, 
                                    k = 6,
                                    bs = "cc",
                                    fx = FALSE,
                                    m = 1,
                                    pc = NULL),
                                   data = fates,
                                   knots = as.list(unname(knot.list)))

smooth.test.1$BD[1, 1] * fates.1$week
basis

test <- predict(ml.model, type = "lpmatrix")


head(-0.05681786 / test[, 2])

test %*% ml.model$coefficients

sweep(test, 2, ml.model$coefficients, `*`)

plot(fates.1$week, exp(rowSums(sweep(test, 2, ml.model$coefficients, `*`))))

# 10 Jan 2025
# I have no idea why the gam() function decreases the basis functions by 1
# seems like overall I'm on the right track, I just need to deal with priors
