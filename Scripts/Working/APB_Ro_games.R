library(tidyverse)
library(readr)

hd <- read_csv("./Data/PopulationGrowth.csv")
glimpse(hd)
filter(hd, Stage =="B0") %>% select(Survival) %>% table()

hd %>% filter(Stage == "B1"|Stage == "B2")

Ro <- hd %>%
  filter(Stage == "B1"|Stage == "B2") %>% 
  mutate(lxmx = Survival * Mean.Rep) %>% 
  ungroup() %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  summarise(Ro = sum(lxmx, na.rm = TRUE))

ggplot(Ro, aes(x = Temperature, y = Ro, group = Treatment, colour = Treatment))+
  geom_point()+
  geom_smooth(span = 2, se = FALSE)+
  facet_grid(Experiment~Clone)

library(lme4)
library(car)
library(MESS)

mod <- lmer(Ro ~ poly(Temperature,2)*Treatment*Experiment + 
              (poly(Temperature,2)+Treatment+Experiment|Clone), data = Ro)
summary(mod)  
Anova(mod)

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(Ro$Treatment),
  Experiment = unique(Ro$Experiment),
  Clone = unique(Ro$Clone)
)

# predicted Ro ----
fixed_pred <- predict(mod, newdata = newX,re.form = NA)
clone_pred <- predict(mod, newdata = newX, re.form = ~(Experiment+Treatment+poly(Temperature,2)|Clone))


# plot data Ro
pd <- data.frame(newX, fixed_pred, clone_pred)

# graph the curves Ro
ggplot(pd, aes(x = Temperature, y = fixed_pred))+
  geom_line(size = 2)+
  # geom_line(aes(x = Temperature, y = clone_pred, colour = Clone), 
  #            size = 1, alpha = 0.5)+
  #geom_jitter(data = Ro, aes(x = Temperature, y = Ro.B2, colour = Clone), alpha = 0.3,
  #            height = 0, width = 1)+
  facet_grid(Experiment ~ Treatment)+
  labs(y = "Population Growth Rate (Ro)")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

# AUC Ro ----
# calcuate area under the curves
AUC <- pd %>% group_by(Clone, Treatment, Experiment) %>% 
  summarise(AUC = auc(x = Temperature, y = clone_pred - min(clone_pred)))

AUCsum <- AUC %>% group_by(Treatment, Experiment) %>% 
  summarise(meanAUC = mean(AUC),
            seAUC = sd(AUC)/sqrt(sum(!is.na(AUC))))

# predation increases marginally the AUC.
ggplot(AUCsum, aes(x = Treatment, y = meanAUC, ymin = meanAUC - seAUC, ymax = meanAUC + seAUC,
                   colour = Experiment, group = Experiment))+
  geom_point(size = 5,position = position_dodge(0.25))+
  geom_line(position = position_dodge(0.25))+
  geom_errorbar(width = 0.1, position = position_dodge(0.25))+
  theme_bw(base_size = 25)

aucMod <- lm(AUC ~ Treatment * Experiment, data = AUC)
Anova(aucMod)

Topts <- pd %>% select(Experiment, Treatment, Temperature, fixed_pred) %>% 
  distinct() %>% 
  group_by(Experiment, Treatment) %>% 
  filter(fixed_pred == max(fixed_pred))

ggplot(pd, aes(x = Temperature, y = fixed_pred, 
               group = Treatment, colour = Treatment))+
  geom_line(size = 2)+
  geom_vline(aes(xintercept = Temperature), data = Topts, col = 'red')+
  geom_hline(aes(yintercept = fixed_pred), data = Topts, col = 'blue')+
  labs(y = "Predicted Population Growth Rate (mm/day)")+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)



# More about SHAPE of the Curves to augment AUC analysis ====

# first get the rows of the predictions
# where the prediction is highest (the peak)
Ro_Max_Topt <- Opt <- pd %>% 
  group_by(Treatment, Experiment) %>% 
  filter(fixed_pred == max(fixed_pred))

# now get the rows at 13C (the min temp)
Ro_TMin <- pd %>% 
  group_by(Treatment, Experiment) %>% 
  filter(Temperature == 13)

# and the rows at 28C (the max temp)
Ro_TMax <- pd %>% 
  group_by(Treatment, Experiment) %>% 
  filter(Temperature == 28)

# put them together and create the Metric column
out <- bind_rows(Ro_Max_Topt, Ro_TMin, Ro_TMax) %>% ungroup() %>% 
  mutate(., Metric = rep(c("Opt", "Min", "Max"), each = 4))

# SHOW how the Shapes change via alterations in values of Ro at Max Temp, Min Temp and Opt Temp.

# first show how values of Ro change among the treatments at Max, Min and Optimal Temp.
TPC1 <- ggplot(out, aes(x = Treatment, y = fixed_pred, group = Experiment, colour = Experiment))+
  geom_point()+geom_line()+
  ylab("Predicted Ro")+
  facet_grid( ~ Metric)+
  theme_bw(base_size = 15)
