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

full_join(tempInd, tempFec)

fec_ind <- full_join(scaleInd, Fec_scale)
glimpse(fec_ind)

filter(fec_ind, Clone == "LD33", Treatment == "Control", Experiment == "Acute", Temperature == 13)

ggplot(fec_ind, aes(x = maxInduction, y = Fec))+
  geom_smooth(method = lm)+
  geom_jitter(aes(colour = Clone))+
  facet_grid(Experiment~Temperature)
