# 20.12.18
# Make All Traits from Raw Data ----

# SGR calculated as linear regression on log Body Size from Intial to Final Juvenile Instar
# Fec calcuated as sum of the number of neonates
# Induction calcuated as max induction
# NO SCALING because 0 matters to TPC.

#libraries
library(tidyverse)

# data
hebe_all <- read_csv("./Data/AllResults_18_07_18.csv")

# Make the trait data ----

# Somatic Growth Data ----
SGR <- hebe_all %>% 
  filter(Clone != "D10_74", Clone != "LD35", Clone != "LD34") %>% 
  mutate(Experiment = case_when(
    Experiment == "Acute" ~ "Acute",
    Experiment == "Aclim" ~ "Acclim",
    Experiment == "Acclim" ~ "Acclim")) %>% 
  filter(Mature == "J") %>% 
  select(ID, Clone, Temperature, Treatment, Experiment, Age, Body) %>% 
  mutate(log.Body = log(Body)) %>% 
  na.omit()

unique(SGR$Clone)

# calculate GR
SGR_scale <- SGR %>%  
  group_by(ID, Clone, Temperature, Treatment, Experiment) %>% 
  nest() %>% 
  mutate(model = map(.x = data , .f = ~lm(log.Body ~ Age, data = .))) %>% 
  mutate(GR = map_dbl(.x = model, .f = ~coef(.)[2]))

# Fecundity Data ----

Fec_scale <- hebe_all %>% 
  filter(Clone != "D10_74", Clone != "LD35", Clone != "LD34", Clone != "D10_74") %>% 
  mutate(Experiment = case_when(
    Experiment == "Acute" ~ "Acute",
    Experiment == "Aclim" ~ "Acclim",
    Experiment == "Acclim" ~ "Acclim")) %>% 
  filter(Mature == "A"&No.Neonates>=1) %>% 
  group_by(ID, Clone, Temperature, Treatment, Experiment) %>% 
  mutate(Fec = sum(No.Neonates, na.rm = TRUE)) %>%
  select(Clone, Temperature, Treatment, Experiment, Fec)

# Induction Data
maxInd_scale <- hebe_all %>%
  filter(Clone != "D10_74", Clone != "LD35", Clone != "LD34") %>% 
  mutate(Experiment = case_when(
    Experiment == "Acute" ~ "Acute",
    Experiment == "Aclim" ~ "Acclim",
    Experiment == "Acclim" ~ "Acclim")) %>% 
  filter(Mature == "J") %>% 
  select(ID, Clone, Temperature, Treatment, Experiment, Spikes, Pedestal) %>% 
  mutate(IndScore = Spikes*10 + Pedestal) %>% 
  group_by(ID, Clone, Temperature, Treatment, Experiment) %>% 
  mutate(maxInd = max(IndScore)) %>% 
  select(Clone, Temperature, Treatment, Experiment, maxInd)


# Initial Plots of three traits by temperature ----
p_Gr <- ggplot(SGR_scale, aes(x = Temperature, y = GR,
                              colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE)+
  facet_grid(Experiment ~ Clone)

p_Fec <- ggplot(Fec_scale, aes(x = Temperature, y = Fec,
                               colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE)+
  facet_grid(Experiment ~ Clone)

p_Ind <- ggplot(maxInd_scale, aes(x = Temperature, y = maxInd,
                                  colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Clone)

gridExtra::grid.arrange(p_Gr, p_Fec, p_Ind)
