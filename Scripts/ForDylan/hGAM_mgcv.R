library(tidyverse)
library(splines)
library(mgcv)

## Source the Trait making script ----

# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are scaled.
# it will produce a plot of all three traits
# takes a few seconds.

source("./Scripts/Working/MakeAllTraits_APB.R")

# look at the data
glimpse(Fec_scale)

## further process the data for gam() - ordered factors etc ----

# un-ordered factors needed and
# temperature needs to be numeric strictly
# create 4 level factor of treatment (2 levels) and experiment (2 levels)
# create 160 level factor of clone (~10), treatment, experiment combos

Fec_scale <- Fec_scale %>% ungroup %>% 
  mutate(Clone = factor(Clone, ordered = FALSE)) %>% 
  mutate(Temperature = as.numeric(Temperature)) %>% 
  # create factor that represents experimental treatments
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), 
                        ordered = FALSE)) %>% 
  # create factor that represents clone and experimental treatments
  # use unite?
  # use "-" for separate later, as "." and "_" are messy
  mutate(ClTrExp = factor(paste(Clone, Treatment,Experiment, sep = "-"), 
                        ordered = FALSE))

## Working with the fectundity data ----
ggplot(Fec_scale, aes(x = Temperature, y = Fec,
                      group = Clone, colour = Clone))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_wrap(~TrExp, ncol = 2)

#plot by 160 level factor: Clone-Treatment-Experiment
ggplot(Fec_scale, aes(x = Temperature, y = Fec, group = ClTrExp))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = F)

## MGCV gam models with Fecundity ----
## see peerJ: https://peerj.com/preprints/27320v1/#supp-1

# models here either have 40 level Clone-Treatment-Experiment Factor
# or 4 level Treatment - Experiment factor with 10 levels of clone

## Model with 40 level Clone, Treatment, Experiment factor ----

# follows first class of model and examples in peerJ
# Model 3 specification: group level trends different smoothers
# A single common smoother plus group-level smoothers
# with DIFFERENT wiggliness
# variation revealed, but not much
# does not use 'experimental design'
# m = 1 creates strong penalty and limited smooth
# setting k = 3 creates smoothest picture... 
# higher k generate limited smoothness and Odd shapes

mod3_ClTrExp <- gam(Fec ~ s(Temperature, bs = "tp", k = 3, m = 2)+
                     s(Temperature, by = ClTrExp, bs = "ts", k = 3, m = 1)+
                     s(ClTrExp, bs = "re", k = 40),
                   data = Fec_scale,
                   drop.unused.levels = FALSE, 
                   method="REML")

# # The simpler model: A single common smoother plus 
# group-level smoothers that have the SAME wiggliness (throws error)
# mod2_ClTrExp <- gam(Fec ~ s(Temperature, k = 3, m = 2)+
#                       s(Temperature, ClTrExp, bs = "fs", k = 3, m = 2),
#                     data = Fec_scale,
#                     drop.unused.levels = FALSE, 
#                     method="REML")
# 
# # more complex justified
# AIC(mod3_ClTrExp) - AIC(mod2_ClTrExp)

# plot(mod_ClTrExp) # too big for laptop screen

## Plotting random effects detail
newX <- expand.grid(Temperature = seq(from = 13, to = 28, by = 0.5),
                    ClTrExp = unique(Fec_scale$ClTrExp))

newY_3 <- predict(mod3_ClTrExp, newdata = newX, type = "response")

plotThese <- data.frame(newX, fittedFec_3 = newY_3) %>% 
  # pop the factor back into separated groups for plotting
  separate(col = ClTrExp, 
           into = c("Clone","Treatment","Experiment"), 
           sep = "-")

# plot model 3 with fully grouped factor
ggplot(plotThese, aes(x = Temperature, y = fittedFec_3,
                      group = Clone, colour = Clone))+
  geom_line()+
  # predictions beyond data for some clones go neg (missing)
  ylim(0,50)+
  facet_grid(Experiment ~ Treatment)

## ALTERNATIVE approach - mapping onto Zooplankton example peerJ ----
# treat clone as random, TrExp as factor
# follows idea of NO common global smoother; 

# The Complex Model (model 5)
# group level trends with DIFFERENT smoothers
# gives different baseline (mean) to each treatment
# and alows differences among clones within treatments
# does not reveal much among clone variation in wiggle
# might be what we want
# not sure what to make m =  and k = 
mod5_TrExp <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 5)+
                s(TrExp, by = Clone, bs = "re")+
                s(Temperature, by = Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

# Model4: group level trends, SAME smoothers
mod4_TrExp <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 5)+
                s(Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

## complex model justified by lower AIC
AIC(mod5_TrExp) - AIC(mod4_TrExp)

# plotting core baselines for each treatment
par(mfrow = c(2,2))
plot(mod5_TrExp, select = 1, scale = 0, ylab = "Control - Acclim",
     shade = TRUE)
plot(mod5_TrExp, select = 2, scale = 0, ylab = "Control - Acute",
     shade = TRUE)
plot(mod5_TrExp, select = 3, scale = 0, ylab = "Predator - Acclim",
     shade = TRUE)
plot(mod5_TrExp, select = 4, scale = 0, ylab = "Predator - Acute",
     shade = TRUE)

## Plotting random effects detail for TrExp models ----
# does not reveal very high levels of clonal variation
newX <- expand.grid(Temperature = seq(from = 13, to = 28, by = 0.5),
                    Clone = levels(Fec_scale$Clone),
                    TrExp = levels(Fec_scale$TrExp))

newY <- predict(mod5_TrExp, newdata = newX, type = "response")
plotThese <- data.frame(newX, fittedFec = newY)

# need to get core smooth onto plots
ggplot(plotThese, aes(x = Temperature, y = fittedFec,
                            group = Clone, colour = Clone))+
  geom_line()+
  facet_wrap(~TrExp, ncol = 2)


