library(tidyverse)
library(splines)
library(mgcv)

# source the Trait making script ----
# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are scaled.
# it will produce a plot of all three traits
# takes a few seconds.
source("./Scripts/Working/MakeAllTraits_APB.R")

glimpse(Fec_scale)

# ordered factors needed and
# temperature needs to be numeric strictly
# create 4 level factor
Fec_scale <- Fec_scale %>% ungroup %>% 
  mutate(Clone = factor(Clone, ordered = FALSE)) %>% 
  mutate(Temperature = as.numeric(Temperature)) %>% 
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), 
                        ordered = FALSE))
  
# Quick View of the data again.
ggplot(Fec_scale, aes(x = Temperature, y = Fec,
                      group = Clone, colour = Clone))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_wrap(~TrExp, ncol = 2)

## MGCV gam models with Fecundity ----

## These seem to be the models to use... complex model is favoured
# This works - I think this might be the good model to use
# Temp x TrExp as baseline
# clone as re
# Group level trends, similar smoothers
mod2_4 <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 4)+
                s(Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_4, pages = 1)

# The More Complex Model
# gives different baseline (mean) to each treatment
# and alows differences among clones within treatments
mod2_5 <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 4)+
                s(TrExp, by = Clone, bs = "re")+
                s(Temperature, by = Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

## complex model justified
anova(mod2_5, mod2_4, test = "Chisq")

# plotting core baselines for each treatment
par(mfrow = c(2,2))
plot(mod2_5, select = 1, scale = 0, ylab = "Control - Acclim",
     shade = TRUE)
plot(mod2_5, select = 2, scale = 0, ylab = "Control - Acute",
     shade = TRUE)
plot(mod2_5, select = 3, scale = 0, ylab = "Predator - Acclim",
     shade = TRUE)
plot(mod2_5, select = 4, scale = 0, ylab = "Predator - Acute",
     shade = TRUE)

## Plotting random effects detail
newX <- expand.grid(Temperature = seq(from = 13, to = 28, by = 0.5),
                    Clone = levels(Fec_scale$Clone),
                    TrExp = levels(Fec_scale$TrExp))

newY <- predict(mod2_5, newdata = newX, type = "response")
plotThese <- data.frame(newX, fittedFec = newY)

# need to get core smooth onto plots
p5 <- ggplot(plotThese, aes(x = Temperature, y = fittedFec,
                            group = Clone, colour = Clone))+
  geom_line()+
  facet_wrap(~TrExp, ncol = 2)

p5


## Plotting model 4 and compare to model 5
newY <- predict(mod2_4, newdata = newX, type = "response")
plotThese <- data.frame(newX, fittedFec = newY)

# need to get core smooth onto plots
p4 <- ggplot(plotThese, aes(x = Temperature, y = fittedFec,
                            group = Clone, colour = Clone))+
  geom_line()+
  facet_wrap(~TrExp, ncol = 2)

gridExtra::grid.arrange(p4, p5, ncol = 2)


# # This is the simplest model, but lacks the Treatments and Experiment
# 
# mod2_1 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
#               s(Temperature, by = Clone, bs = "re"),
#             data = Fec_scale,
#             drop.unused.levels = FALSE, 
#             method="REML")
# 
# plot(mod2_1, pages = 1)
# 
# 
# # This works but....
# # tp spline for base
# # Treatment_Experiment as re
# # Clone nested in Treatment_Experiment as re
# mod2_2 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
#                 s(TrExp, bs = "re")+
#                 s(TrExp, Clone, bs = "re"),
#               data = Fec_scale,
#               drop.unused.levels = FALSE, 
#               method="REML")
# 
# plot(mod2_2, pages = 1)
# 
# # newXs
# newX <- expand.grid(Temperature = seq(13, 28, 0.1),
#                     TrExp = unique(Fec_scale$TrExp),
#                     Clone = unique(Fec_scale$Clone))
# 
# # predictions
# pred_mod2_2 <- predict(mod2_2, se.fit=TRUE, newdata = newX)
# 
# # housekeeping
# plotThese <- transform(newX, fits = pred_mod2_2$fit, fit_se = pred_mod2_2$se.fit)
# 
# # plot
# ggplot(data=plotThese, aes(x=Temperature, y=fits, group=Clone, colour = Clone)) +
#   facet_wrap(~TrExp) +
#   geom_line() +
#   geom_point(data = Fec_scale, aes(x = Temperature, y = Fec, 
#                                    group = Clone, colour = Clone))+
#   geom_ribbon(aes(ymin=fits-2*fit_se,
#                   ymax=fits+2*fit_se, fill = Clone), alpha=0.1, colour = NA)
# 
# 
# # this doesn't work
# # tp for base
# # Temperature by TrExp as factor spline
# # Temperature by clone as re
# mod2_3 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
#               s(Temperature, by = TrExp, bs = "fs")+
#               s(Temperature, by = Clone, bs = "re"),
#             data = Fec_scale,
#             drop.unused.levels = FALSE, 
#             method="REML")
# 
# 
# 
# 
#   
