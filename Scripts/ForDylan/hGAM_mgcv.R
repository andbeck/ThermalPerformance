library(tidyverse)
library(splines)
library(mgcv)

# data
hd <- read_csv("./Data/PopulationGrowth.csv")
glimpse(hd)


# Get Ro data in order
Ro <- hd %>%
  filter(Stage == "B1"|Stage == "B2") %>% 
  mutate(lxmx = Survival * Mean.Rep) %>% 
  ungroup() %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  summarise(Ro = sum(lxmx, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # create combined factor, unordered
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), ordered = FALSE)) %>% 
  # ensure clone is unordered factor (from paper)
  mutate(Clone = factor(Clone, ordered = FALSE)) %>% 
  # ensure Temperature is numeric
  mutate(Temperature = as.numeric(Temperature))

# first view of it
ggplot(Ro, aes(x = Temperature, y = Ro, group = Clone, colour = Clone))+
  geom_point()+
  geom_smooth(span = 2, se=FALSE)+
  facet_grid(Treatment ~ Experiment)


# MGCV gam models
# This works but lacks the Treatments and Experiment.
mod2_1 <- gam(Ro ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = Clone, bs = "re"),
            data = Ro,
            drop.unused.levels = FALSE, 
            method="REML")
plot(mod2_1, pages = 1)


# This works but....
# tp spline for base
# Treatment_Experiment as re
# Clone nested in Treatment_Experiment as re
mod2_2 <- gam(Ro ~ s(Temperature, bs = "tp", k = 4)+
                s(TrExp, bs = "re")+
                s(TrExp, Clone, bs = "re"),
              data = Ro,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_2, pages = 1)

# this doesn't work
# tp for base
# Temperature by TrExp as factor spline
# Temperature by clone as re
mod2_3 <- gam(Ro ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = TrExp, bs = "fs")+
              s(Temperature, by = Clone, bs = "re"),
            data = Ro,
            drop.unused.levels = FALSE, 
            method="REML")


# This works
# Temp x TrExp as baseline
# clone as re
mod2_4 <- gam(Ro ~ s(Temperature, by = TrExp, bs = "tp", k = 4)+
                s(Clone, bs = "re"),
              data = Ro,
              drop.unused.levels = FALSE, 
              method="REML")

plot(mod2_4, pages = 1)

# This does not work either, matches format at end of PeerJ
# Temp for base
# Temp by TrExo as factor spline
# TrExp with Clone nested as re
mod2_5 <- gam(Ro ~ s(Temperature, bs = "tp", k = 4)+
              s(Temperature, by = TrExp, bs = "fs")+
              s(TrExp, by = Clone, bs = "re"),
            data = Ro,
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
