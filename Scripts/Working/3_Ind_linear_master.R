library(tidyverse)
library(lme4)
library(car)
library(heplots)
library(MESS)

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")

# Create a scaled dataset and a response one with Factorial ----------------------------
# Omit D10_A13 because so much missing

scaleDat <- SGR_Ind %>% as.data.frame() %>% 
  mutate(maxInduction = scale(maxInduction),
         factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
  na.omit() %>% 
  data.frame()

scaleDat %>% 
  group_by(Temperature,Experiment, Treatment, Clone) %>% 
  summarise(count = sum(!is.na(maxInduction))) %>% 
  filter(count == 3)

# first plot of scaled data
ggplot(scaleDat, aes(x = Temperature, y = maxInduction, 
                     colour = Treatment))+
  geom_jitter()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  ylab("Induction Score (Scaled)")+
  theme_bw()

ggplot(scaleDat, aes(x = Temperature, y = maxInduction, group = Experiment))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE, size = 3, colour = 'black')+
  geom_smooth(aes(group = factorial), method = lm,
              se = FALSE, size = 1, colour = 'grey')+
  facet_wrap(~Experiment)+
  theme_bw()

# Model --------------------------------
# FULL random effects structure justified by LRT
mod <- lmer(maxInduction ~ Temperature * Treatment * Experiment +
              (Temperature|Clone),
            data = scaleDat)

# mod <- lmer(maxInduction ~ Temperature * Treatment * Experiment +
#               (Temperature+Treatment+Experiment|Clone),
#             data = scaleDat)

# anova(mod, mod2)

summary(mod)
Anova(mod, test.statistic = "F")

# PREDICTIONS

# new X's grid
newX <- expand.grid(
  Temperature = seq(13,28,length = 500),
  Treatment = unique(scaleDat$Treatment),
  Experiment = unique(scaleDat$Experiment),
  Clone = unique(scaleDat$Clone)
)

# predicted Fecundities
fixed_pred <- predict(mod, newdata = newX,re.form = NA)
clone_pred <- predict(mod, newdata = newX, re.form = ~Temperature|Clone)

# housekeeping
pd <- data.frame(newX, fixed_pred, clone_pred) %>% 
  mutate(Experiment = factor(Experiment, levels = c("Acclim", "Acute")))

# graph the curves Fecundity
ggplot(pd, aes(x = Temperature, y = fixed_pred))+
  geom_line(size = 2)+
  geom_line(aes(x = Temperature, y = clone_pred, colour = Clone), 
            size = 1, alpha = 0.5)+
  # geom_jitter(data = Ro, aes(x = Temperature, y = Ro.B2, colour = Clone), alpha = 0.3,
  #             height = 0, width = 1)+
  facet_grid(Experiment ~ Treatment)+
  labs(y ="Induction Score (Scaled)")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

