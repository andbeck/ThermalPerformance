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

# Fecundity Data ----

Fec_scale <- hebe_all %>% 
  filter(Clone != "D10_74", Clone != "LD35", Clone != "LD34") %>% 
  mutate(Experiment = case_when(
    Experiment == "Acute" ~ "Acute",
    Experiment == "Aclim" ~ "Acclim",
    Experiment == "Acclim" ~ "Acclim")) %>% 
  filter(Mature == "A"&No.Neonates>=1) %>% 
  group_by(ID, Clone, Temperature, Treatment, Experiment) %>% 
  mutate(Fec = sum(No.Neonates, na.rm = TRUE)) %>%
  select(Clone, Temperature, Treatment, Experiment, Fec) %>% 
  distinct()

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")

# Create a scaled dataset and a response one with Factorial ----------------------------
# Omit D10_A13 because so much missing

SGR_scale <- SGR_Ind %>% as.data.frame() %>% 
  filter(Clone != "D10_A13") %>% 
  mutate(maxInduction = scale(maxInduction),
         Growth0 = scale(Growth0),
         Growth1 = scale(Growth1),
         factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
  na.omit() %>% 
  data.frame()

# scaleDat %>% 
#   group_by(Temperature,Experiment, Treatment, Clone) %>% 
#   summarise(count = sum(!is.na(Growth0))) %>% 
#   filter(count == 1)

maxInd_scale <- SGR_Ind %>% as.data.frame() %>% 
  mutate(maxInd = scale(maxInduction),
         factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
  na.omit() %>% 
  data.frame()


# Initial Plots of three traits by temperature ----
p_Gr <- ggplot(SGR_scale, aes(x = Temperature, y = Growth0,
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
