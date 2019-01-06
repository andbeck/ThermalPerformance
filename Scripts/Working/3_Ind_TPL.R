library(tidyverse)
library(splines)
library(lme4)
#devtools::install_github('ZheyuanLi/SplinesUtils')
library(SplinesUtils)
library(broom)
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
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Clone)+
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
              (Temperature+Treatment+Experiment|Clone),
            data = scaleDat)

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

# predictions at level 1 (all Random Effects)
pred_full <- predict(mod, newdata = newX)
pred_fix <- predict(mod, newdata = newX, re.form = NA)

# housekeeping
predDat <- data.frame(newX, pred_full, pred_fix)

# BUILD AND VISUALISE  EFFECTS VIEW ====================================

# visualise FIXED responses
ggplot(predDat, aes(x = Temperature, y = pred_fix,
                    group = Treatment, colour = Treatment))+
  geom_line(size = 2)+
  labs(y = "max Induction Score")+
  scale_colour_manual(values = c(Control = "Black", Predator = "Red"))+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)

# with RE
ggplot(predDat, aes(x = Temperature, y = pred_full, 
                    group = Clone, colour = Clone))+
  geom_line()+
  geom_line(aes(x = Temperature, y = pred_fix), colour = "grey20")+
  facet_grid(Treatment ~ Experiment)
