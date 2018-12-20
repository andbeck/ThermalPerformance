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
Fec_scale <- Fec_scale %>% 
  mutate(Clone = factor(Clone, ordered = FALSE)) %>% 
  mutate(Temperature = as.numeric(Temperature)) %>% 
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), 
                        ordered = FALSE))
  
## MGCV gam models with Fecundity ----

# This works but lacks the Treatments and Experiment

mod2_1 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = Clone, bs = "re"),
            data = Fec_scale,
            drop.unused.levels = FALSE, 
            method="REML")

plot(mod2_1, pages = 1)


# This works but....
# tp spline for base
# Treatment_Experiment as re
# Clone nested in Treatment_Experiment as re
mod2_2 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
                s(TrExp, bs = "re")+
                s(TrExp, Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_2, pages = 1)

# this doesn't work
# tp for base
# Temperature by TrExp as factor spline
# Temperature by clone as re
mod2_3 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = TrExp, bs = "fs")+
              s(Temperature, by = Clone, bs = "re"),
            data = Fec_scale,
            drop.unused.levels = FALSE, 
            method="REML")


# This works
# Temp x TrExp as baseline
# clone as re
mod2_4 <- gam(Fec ~ s(Temperature, by = TrExp, bs = "tp", k = 4)+
                s(Clone, bs = "re"),
              data = Fec_scale,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_4, pages = 1)

# This does not work either, matches format at end of PeerJ
# Temp for base
# Temp by TrExo as factor spline
# TrExp with Clone nested as re
mod2_5 <- gam(Fec ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = TrExp, bs = "fs")+
              s(TrExp, by = Clone, bs = "re"),
            data = Fec_scale,
            drop.unused.levels = FALSE, 
            method="REML")



# prediction plots
# need to get core smooth as well
newX <- expand.grid(Temperature = seq(from = 13, to = 28, by = 0.5),
                    Clone = levels(Ro$Clone),
                    TrExp = levels(Ro$TrExp))

newY <- predict(mod2_4, newdata = newX, type = "response")
plotThese <- data.frame(newX, fittedRo = newY)

# need to get core smooth onto plots
ggplot(plotThese, aes(x = Temperature, y = fittedRo,
                      group = Clone, colour = Clone))+
  geom_line()+
  facet_wrap(~TrExp, ncol = 2)
  
plot(mod2_4) 
