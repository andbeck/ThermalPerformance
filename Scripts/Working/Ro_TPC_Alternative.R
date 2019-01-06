library(tidyverse)
library(lme4)
library(car)
library(heplots)
library(MESS)

# data
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

## Ro ------------------------------------------

# Create a scaled dataset and a response one with Factorial

scaleDat <- Ro %>% as.data.frame() %>% 
  mutate(Ro = scale(Ro),
         factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
  na.omit() %>% 
  data.frame()

# first plot of scaled data
ggplot(scaleDat, aes(x = Temperature, y = Ro, 
                     colour = Treatment))+
  geom_jitter()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  ylab(expression(paste("Reproduction (", Sigma, "lxmx)")))+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  theme_bw()

#Acute Screws up reproduction!
ggplot(scaleDat, aes(x = Temperature, y = Ro, group = Experiment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, size = 3, colour = 'black')+
  geom_smooth(aes(group = factorial), method = lm, 
              formula = y ~ poly(x, 2), se = FALSE, size = 1, colour = 'grey')+
  facet_wrap(~Experiment)+
  theme_bw()

## bs method VERSION XFILE ---------------------------------------

# fit model using scaled data and complex RE
# or simple RE

# experiment|Clone not justified by LRT
# temperature|Clone NOT justified by LRT
# bs(Temperature)|Clone justified
# bs(Temeprature) + Experiment NO
# bs(Temeprature) + Experiment + Treatment No
# library(splines)
# mod <- lmer(Ro ~ bs(Temperature) * Treatment * Experiment +
#                (bs(Temperature)|Clone),
#              data = scaleDat)

# NOTE that bs() function for Acute group creates VERY awkward shape
# poly(T,3) produces same as bs() - dafulting to knots 3, clearly
# result is the same though (qual) as the peack performance is at v. low temps
mod <- lmer(Ro ~ poly(Temperature, 2) * Treatment * Experiment +
              (1|Clone),
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
  labs(y = "Predicted Ro")+
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
