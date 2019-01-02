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
  # create factor that represents experimental treatments
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), 
                        ordered = FALSE)) %>% 
  # create factor that represents clone and experimental treatments
  # use unite?
  mutate(ClTrExp = factor(paste(Clone, Treatment,Experiment, sep = "-"), 
                        ordered = FALSE))

# Quick View of the data again.
ggplot(Fec_scale, aes(x = Temperature, y = Fec,
                      group = Clone, colour = Clone))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_wrap(~TrExp, ncol = 2)

#plot by 'collapsed factor'
ggplot(Fec_scale, aes(x = Temperature, y = Fec, group = ClTrExp))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = F)

## MGCV gam models with Fecundity ----

## These are one version of the models to use... 
## treat clone as random, TrExp as factor

# The More Complex Model (model 5)
# group level trends with different smoothers
# gives different baseline (mean) to each treatment
# and alows differences among clones within treatments
# does not reveal much among clone variation in wiggle
# might be what we want
# not sure what to make m =  and k = 
mod2_5 <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 5)+
                s(TrExp, by = Clone, bs = "re")+
                s(Temperature, by = Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

# Model4: group level trends, similar smoothers
mod2_4 <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 4)+
                s(Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_4, pages = 1)

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
# does not reveal very high levels of clonal variation
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

## Model with Clone, Treatment, Experiment grouping variable
# Model 3 specification here: group level trends different smoothers
# but single baseline
# higher variation revealed, but not much
# does not use 'experimental design'
# m = 1 creates strong penalty and limited smooth
# setting k = 3 creates smoothest picture... 
# higher k generate limited smoothness and Odd shapes
mod_ClTrExp <- gam(Fec ~ s(Temperature, bs = "tp", k = 3, m = 2)+
                     s(Temperature, by = ClTrExp, bs = "ts", k = 3, m = 1)+
                     s(ClTrExp, bs = "re", k = 40),
                     data = Fec_scale,
                     drop.unused.levels = FALSE, 
                     method="REML")


plot(mod_ClTrExp)

## Plotting random effects detail
newX <- expand.grid(Temperature = seq(from = 13, to = 28, by = 0.5),
                    ClTrExp = unique(Fec_scale$ClTrExp))

newY <- predict(mod_ClTrExp, newdata = newX, type = "response")
plotThese <- data.frame(newX, fittedFec = newY) %>% 
  separate(col = ClTrExp, into = c("Clone","Treatment","Experiment"), sep = "-")

#need to get core smooth onto plots
ggplot(plotThese, aes(x = Temperature, y = fittedFec,
                            group = Clone, colour = Clone))+
  geom_line()+
  ylim(0,50)+
  facet_grid(Experiment ~ Treatment)
