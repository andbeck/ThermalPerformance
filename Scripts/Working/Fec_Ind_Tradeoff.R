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
  na.omit() %>% 
  select(Clone, Temperature, Treatment, Experiment, maxInduction)

# Fec Data
source("./Scripts/Working/MakeFecData_APB.R")
glimpse(Fec_scale)

Fec_scale <- Fec_scale %>% 
  select(Clone, Temperature, Treatment, Experiment, Fec)


# look at the data
glimpse(scaleInd)
filter(scaleInd, Clone == "LD33", Treatment == "Control", Experiment == "Acute", Temperature == 13) %>% data.frame()
glimpse(Fec_scale)
filter(Fec_scale, Clone == "LD33", Treatment == "Control", Experiment == "Acute", Temperature == 13) %>% data.frame()

fec_ind <- left_join(scaleInd, Fec_scale)
glimpse(fec_ind)
filter(fec_ind, Clone == "LD33", Treatment == "Control", Experiment == "Acute", Temperature == 13)

ggplot(fec_ind, aes(x = maxInduction, y = Fec))+
  geom_smooth(method = lm)+
  geom_jitter(aes(colour = Clone))+
  facet_grid(Experiment~Temperature)
