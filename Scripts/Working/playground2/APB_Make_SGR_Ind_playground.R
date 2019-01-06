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
  na.omit() %>% 
  distinct()

unique(SGR$Clone)

# calculate GR using slope or regression of log.Body vs. Age
# not what Hebe did
SGR_scale <- SGR %>%  
  group_by(ID, Clone, Temperature, Treatment, Experiment) %>% 
  nest() %>% 
  mutate(model = map(.x = data , .f = ~lm(log.Body ~ Age, data = .))) %>% 
  mutate(GR = map_dbl(.x = model, .f = ~coef(.)[2]))


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
  select(Clone, Temperature, Treatment, Experiment, maxInd) %>% 
  distinct()
