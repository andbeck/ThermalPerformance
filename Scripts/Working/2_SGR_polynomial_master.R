# Somatic Growth Rate Analysis
# uses Hebe's data and Growth0 (size at maturity) endpoint

# Take home messages

# Significant effects of Treatment and Experiment
# acclimation matters and predation matters
# No interactions justified

# effect size of AUC INTERACTION moderate 
# Suggests Fast-Slow tied to acclimation and predation

# Very small effect sizes of Topt and Popt
# NOT generalist specialist or hot cold.

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
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  ylab("Somatic Growth Rate (mm/day)")+
  theme_bw()

ggplot(scaleDat, aes(x = Temperature, y = Growth0, group = Experiment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, size = 3, colour = 'black')+
  geom_smooth(aes(group = factorial), method = lm, 
              formula = y ~ poly(x, 2), se = FALSE, size = 1, colour = 'grey')+
  facet_wrap(~Experiment)+
  theme_bw()

## SGR lmer ---------------------------------------

mod <- lmer(Growth0 ~ poly(Temperature,2) * Treatment * Experiment +
              (poly(Temperature,2)|Clone),
            data = scaleDat)

summary(mod)
Anova(mod, test.statistic = "F")

# Plot the results ----

# showing variation among clones
# and fixed effect

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(scaleDat$Treatment),
  Experiment = unique(scaleDat$Experiment),
  Clone = unique(scaleDat$Clone)
)

# predicted Fecundities
fixed_pred <- predict(mod, newdata = newX,re.form = NA)
clone_pred <- predict(mod, newdata = newX, 
                      re.form = ~(poly(Temperature,2)|Clone))

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
  labs(y = "Somatic Growth Rate (mm/day)")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

## AUC and Popt/Topt analysis ----

# calcuate area under the curves ----
AUC <- pd %>% group_by(Clone, Treatment, Experiment) %>%
  summarise(AUC = auc(x = Temperature, y = clone_pred - min(clone_pred)))

AUCsum <- AUC %>% group_by(Treatment, Experiment) %>% 
  summarise(meanAUC = mean(AUC),
            seAUC = sd(AUC)/sqrt(sum(!is.na(AUC))))

# GRAPH AUC: predation increases the AUC ----
ggplot(AUCsum, aes(x = Experiment, y = meanAUC, ymin = meanAUC - seAUC, ymax = meanAUC + seAUC,
                   colour = Treatment, group = Treatment))+
  geom_point(size = 5,position = position_dodge(0.25))+
  geom_line(position = position_dodge(0.25))+
  geom_errorbar(width = 0.1, position = position_dodge(0.25))+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  ylab("AUC Somatic Growth")+
  theme_bw(base_size = 15)

# MODEL AUC ----
aucMod <- lm(AUC ~ Treatment * Experiment, data = AUC)
Anova(aucMod)
heplots::etasq(aucMod, anova = TRUE) #substantial effect size on interaction

# T- and P-opt values for analysis ----
# use clone_pred

P_T_opts <- pd %>% group_by(Clone, Treatment, Experiment) %>% 
  # get the rows where clone preds are max in each group
  filter(clone_pred == max(clone_pred)) %>% 
  select(Temperature, Treatment, Experiment, Clone, clone_pred) %>% 
  # rename to conventions
  rename(Topt = Temperature, Popt = clone_pred)

# Topt/Popt analysis ----
mod_Topts <- lm(Topt ~ Treatment * Experiment, data = P_T_opts)
mod_Popts <- lm(Popt ~ Treatment * Experiment, data = P_T_opts)
heplots::etasq(mod_Topts, anova = TRUE) # small effect sizes on Topt
heplots::etasq(mod_Popts, anova = TRUE) # small effect sizes Performance

# T- and P-opt values for plotting (use fixed_pred) ----

graph_PT <- pd %>% ungroup() %>% 
  select(Temperature, Treatment, Experiment, fixed_pred) %>% 
  group_by(Treatment, Experiment) %>% 
  filter(fixed_pred == max(fixed_pred)) %>%
  rename(Topt = Temperature, Popt = fixed_pred) %>% 
  distinct()

# GRAPH Topt/Popt ----
# a thing of beauty

ggplot(pd, aes(x = Temperature, y = fixed_pred, 
               group = Treatment, colour = Treatment))+
  # add the background clone lines
  geom_vline(aes(xintercept = Topt, colour = Treatment), data = P_T_opts, 
             alpha = 0.1)+
  geom_hline(aes(yintercept = Popt, colour = Treatment), data = P_T_opts,
             alpha = 0.1)+
  # add the fixed lines
  geom_vline(aes(xintercept = Topt, colour = Treatment), data = graph_PT, 
             alpha = 0.6)+
  geom_hline(aes(yintercept = Popt, colour = Treatment), data = graph_PT,
             alpha = 0.6)+
  # add the curves on top of everything
  geom_line(size = 1)+
  # the rest
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  labs(y = "Somatic Growth Rate (mm/day)")+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)
