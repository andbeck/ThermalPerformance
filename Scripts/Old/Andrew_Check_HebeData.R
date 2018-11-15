# Hebe All Data
rm(list = ls())

library(tidyverse)

# all data
allDat <- read_csv("AllResults_17_07_18.csv")
allDat$Clone
unique(allDat$Stage)

# check your data and levels of the factors
distinct(allDat, Clone)
distinct(allDat, Treatment)
distinct(allDat, Temperature)
distinct(allDat, Experiment)
distinct(allDat, Stage) # problems
distinct(allDat, Mature) # problems

# fix the factor levels
allDat2 <- 
  allDat %>%
  mutate(Experiment = ifelse(Experiment == "Aclim", "Acclim", Experiment)) %>% 
  mutate(Stage = ifelse(Stage == "Mature", "mt", Stage)) %>% 
  mutate(Clone = ifelse(Clone == "LD34"|Clone == "LD35", "LD33", Clone))

# check
distinct(allDat2, Clone)
distinct(allDat2, Stage)

# get the mature data
matDat <- filter(allDat2, Stage == "mt")

# plot it to isolate other issues
ggplot(matDat, aes(x = factor(Temperature), y = Body, 
                   colour = Clone, shape = Treatment))+
  geom_jitter(height = 0, width = 0.25)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()

# consider removing D10_74
# D10_45 has some very low numbers at 13C


# check sample sizes by temp and treatment
# something seems fishy here.
sampSize <- matDat %>% 
  group_by(Experiment, Treatment, Temperature, Clone) %>% 
  #summarise(nums = sum(!is.na(Body)))
  summarise(nums = n())

ggplot(sampSize, aes(x = factor(Temperature), y = nums, 
                     colour = Treatment, group = Treatment))+
  geom_point()+
  geom_line()+
  facet_grid(Experiment ~ Clone)+
  scale_y_continuous(breaks = 1:10)+
  theme_bw()
