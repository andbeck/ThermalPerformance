library(tidyverse)
library(lme4)
library(car)
library(heplots)
library(MESS)

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")
names(SGR_Ind)

# Create a scaled dataset and a response one with Factorial ----------------------------
# Omit D10_A13 because so much missing

scaleInd <- SGR_Ind %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  mutate(Rep = row_number()) %>% 
  mutate(ID = paste(Rep, Clone, Temperature, Treatment, Experiment, sep = "-")) %>% 
  select(ID, Clone, Temperature, Treatment, Experiment, maxInduction)

# Fec Data
source("./Scripts/Working/MakeFecData_APB.R")
glimpse(Fec_scale)

Fec_scale <- Fec_scale %>% ungroup() %>% 
  select(Clone, Temperature, Treatment, Experiment, Fec) %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  mutate(Rep = row_number()) %>% 
  mutate(ID = paste(Rep, Clone, Temperature, Treatment, Experiment, sep = "-")) %>% 
  select(ID, Clone, Temperature, Treatment, Experiment, Fec)



# look at the data
glimpse(scaleInd)
tempInd <- filter(scaleInd, 
       Clone == "LD33", 
       Treatment == "Control", 
       Experiment == "Acute", 
       Temperature == 13) %>% data.frame()

glimpse(Fec_scale)
tempFec <- filter(Fec_scale, 
       Clone == "LD33", 
       Treatment == "Control", 
       Experiment == "Acute", 
       Temperature == 13) %>% data.frame()

tempFec
tempInd

# CHECK THIS AGAIN.
fec_ind <- full_join(scaleInd, Fec_scale)
glimpse(fec_ind)

# test
filter(fec_ind, 
       Clone == "LD33", Treatment == "Control", 
       Experiment == "Acute", Temperature == 13)

# plot all of the data
ggplot(fec_ind, aes(x = maxInduction, y = Fec))+
  geom_smooth(method = lm, se = FALSE)+
  geom_jitter(aes(colour = Clone))+
  facet_grid(Experiment~Temperature)

# take clone means
mean_fec_ind <- fec_ind %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  summarise(meanInd = mean(maxInduction, na.rm = TRUE),
            meanFec = mean(Fec, na.rm = TRUE))

# factorise by temperature
F0 <-  mean_fec_ind %>% 
  ungroup() %>% 
  mutate(Temperature = factor(Temperature)) %>% 
  filter(Treatment == "Predator")

# plot it
# perhaps a strong ind x temp in acclim, but not in acute?
ggplot(F0, aes(x = meanInd, y = meanFec))+
  # coloured points by temperature
  geom_point(aes(group = Temperature, colour = Temperature))+
  # temperature regressions among clones
  geom_smooth(aes(group = Temperature, colour = Temperature), 
              method = lm, se=FALSE, alpha = 0.5)+
  # overall regression
  geom_smooth(method = lm, se = FALSE, col = 'black')+
  # customise colour scheme
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(Experiment ~ .)

# Models of clone means
modTO <- lm(meanFec ~ meanInd * Temperature * Experiment, data = F0)

Anova(modTO)
summary(modTO)

newX <- expand.grid(
  meanInd = seq(from = 0, to = 100, by = 10),
  Temperature = unique(G0$Temperature),
  Experiment =  unique(G0$Experiment))

newY <- predict(modTO, newdata = newX, interval = 'confidence')

plotThese <- data.frame(newX, newY) %>% 
  rename(predictedSGR = fit, Induction = meanInd)

ggplot(plotThese, aes(x = Induction, y = predictedSGR, 
                      group = Temperature, colour = Temperature))+
  geom_line()+
  geom_point(data = G0, 
             aes(x = meanInd, y = meanSGR))+
  facet_wrap(~Experiment)