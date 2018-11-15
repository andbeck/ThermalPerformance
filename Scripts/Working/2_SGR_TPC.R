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
  filter(Clone != "D10_A13") %>% 
  mutate(maxInduction = scale(maxInduction),
         Growth0 = scale(Growth0),
         Growth1 = scale(Growth1),
         factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
  na.omit() %>% 
  data.frame()

scaleDat %>% 
  group_by(Temperature,Experiment, Treatment, Clone) %>% 
  summarise(count = sum(!is.na(Growth0))) %>% 
  filter(count == 1)

# first plot of scaled data
ggplot(scaleDat, aes(x = Temperature, y = Growth0, 
                     colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()

ggplot(scaleDat, aes(x = Temperature, y = Growth0, group = Experiment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, size = 3, colour = 'black')+
  geom_smooth(aes(group = factorial), method = lm, 
              formula = y ~ poly(x, 2), se = FALSE, size = 1, colour = 'grey')+
  facet_wrap(~Experiment)+
  theme_bw()

## bs method VERSION XFILE ---------------------------------------

# fit model using scaled data and complex RE
# or simple RE

# simple.  Note that (bs(Temperature)|Clone) does not work
# this one has all clones with tunring point
# very similar Topts and variable Popts
mod <- lmer(Growth0 ~ bs(Temperature) * Treatment * Experiment +
              (1|Clone),
            data = scaleDat)

# experiment|Clone justified by LRT
# temperature|Clone NOT justified by LRT
# bs(Temperature)|Clone does not fit
# bs(Temeprature) + Experiment GOOD
# bs(Temeprature) + Experiment + Treatment not better than just with experiment

# RE: this has more variation in Topt.  Some clones
# have no turning point
mod <- lmer(Growth0 ~ bs(Temperature) * Treatment * Experiment +
              (bs(Temperature)+Experiment|Clone),
            data = scaleDat)

# quite big differences with poly(T,2), but not (T,3)
mod <- lmer(Growth0 ~ poly(Temperature,3) * Treatment * Experiment +
              (poly(Temperature,3)+Experiment|Clone),
            data = scaleDat)


summary(mod)
Anova(mod, test.statistic = "F")

# genreate predictions -------------------------------------

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

# BUILD AND VISUALISE FIXED EFFECTS VIEWS of Topt, Popt and AUC ====================================

# Calculate Topts via max of curves
Topts <- predDat %>% select(Experiment, Treatment, Temperature, pred_fix) %>% 
  distinct() %>% 
  group_by(Experiment, Treatment) %>% 
  filter(pred_fix == max(pred_fix))

# Calculate AUCs from curves
AUCs <- predDat %>% select(Experiment, Treatment, Temperature, pred_fix) %>% 
  group_by(Experiment, Treatment) %>%
  summarise(AUC = MESS::auc(x = Temperature, y = pred_fix - min(pred_fix)))
  
# calucate Changes from Topts and Popts
deltas <- left_join(Topts, AUCs) %>% 
  group_by(Experiment) %>% 
  # Control - Predator: neg means C < P
  summarise(deltaT = diff(Temperature),
            deltaP = diff(pred_fix),
            deltaAUC = diff(AUC))

# visualise FIXED responses
ggplot(predDat, aes(x = Temperature, y = pred_fix, group = Treatment, colour = Treatment))+
  geom_line(size = 2)+
  geom_vline(aes(xintercept = Temperature), data = Topts, col = 'grey20', linetype = 'dashed')+
  geom_hline(aes(yintercept = pred_fix), data = Topts, col = 'blue', linetype = 'dashed')+
  labs(y = "Predicted Somatic Growth Rate (mm/day)")+
  scale_colour_manual(values = c(Control = "Black", Predator = "Red"))+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)

# BUILD AND VISUALISE Randome EFFECTS VIEWS of Topt, Popt and AUC ====================================

Topts_clone <- predDat %>% select(Experiment, Treatment, Temperature, Clone, pred_full) %>% 
  distinct() %>% 
  group_by(Experiment, Treatment, Clone) %>% 
  filter(pred_full == max(pred_full))

AUCs_clone <- predDat %>% select(Experiment, Treatment, Clone, Temperature, pred_full) %>% 
  group_by(Experiment, Treatment, Clone) %>%
  summarise(AUC = MESS::auc(x = Temperature, y = pred_full - min(pred_full)))

fullCloneStats <- left_join(Topts_clone, AUCs_clone)
modMan <- lm(cbind(Temperature, pred_full, AUC) ~ Experiment * Treatment, data = fullCloneStats)
Anova(modMan)
summary(modMan)

modTopt <- lm(Temperature ~ Experiment * Treatment, data = fullCloneStats)
etasq(modTopt, anova = TRUE)

modPopt <- lm(pred_full ~ Experiment * Treatment, data = fullCloneStats)
etasq(modPopt, anova = TRUE)

modAUC <- lm(AUC ~ Experiment * Treatment, data = fullCloneStats)
etasq(modAUC, anova = TRUE)

prettyFix <- predDat %>% select(Experiment, Treatment, Temperature, Clone, pred_fix) %>% 
  distinct()

ggplot(predDat, aes(x = Temperature, y = pred_full, group = Clone, colour = Clone))+
  geom_line()+
  geom_line(aes(x = Temperature, y = pred_fix), colour = "grey20")+
  facet_grid(Treatment ~ Experiment)


