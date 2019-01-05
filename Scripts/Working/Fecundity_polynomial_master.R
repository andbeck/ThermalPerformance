# FECUNDITY
# sum of 3-clutches
# Analysis of TPC, AUC and Topt/Popt

# Take home: significant Genetic Variation
# Acclimation vs Acute is substantial
# curves really only under Acclimation
# AUC indicates Fast - Slow
# Topt/Popt indicates Generalist Specialist.


library(tidyverse)
library(splines)
library(lme4)
library(car)
library(MESS)

## Source the Trait making script ----

# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are scaled.
# it will produce a plot of all three traits
# takes a few seconds.

source("./Scripts/Working/MakeAllTraits_APB.R")

## Fecundity Analyses ----

# look at the data
glimpse(Fec_scale)
p_Fec # the plot

# The model ----

# Effect of Temp varies by treatment and by Experiment
# Effect of Treatment does not vary by Experiment (Predation effect is invariant)
mod <- lmer(Fec ~ (poly(Temperature,2)+Treatment+Experiment)^2 + 
              (poly(Temperature,2)|Clone), 
            data = Fec_scale)

summary(mod)  
Anova(mod)

# Plot the results ----

# showing variation among clones
# and fixed effect

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(Fec_scale$Treatment),
  Experiment = unique(Fec_scale$Experiment),
  Clone = unique(Fec_scale$Clone)
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
  labs(y =expression(paste("Fecundity ", Sigma, "3-clutches")))+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

## AUC and Popt/Topt analysis ----

# clones D10_81 and D86_A have non-parabolic and increasing curves
# in control, acute combination (see graph above)
# get rid of these two for the AUC and P_T_opt analyses
# they do not change the qualitative outcome, but generate
# peak performance at highest temperature (possibly true)

pd %>% group_by(Clone, Treatment, Experiment) %>%
  top_n(n=1) %>% filter(Temperature == 28)

pd <- pd %>% group_by(Clone, Treatment, Experiment) %>%
  mutate(ID = paste(Clone, Treatment, Experiment, sep = "_")) %>% 
  filter(ID != "D10_81_Control_Acute", ID != "D86_A_Control_Acute")

# calcuate area under the curves ----
AUC <- pd %>% 
  summarise(AUC = auc(x = Temperature, y = clone_pred - min(clone_pred)))

AUCsum <- AUC %>% group_by(Treatment, Experiment) %>% 
  summarise(meanAUC = mean(AUC),
            seAUC = sd(AUC)/sqrt(sum(!is.na(AUC))))

# GRAPH AUC: predation increases the AUC ----
ggplot(AUCsum, aes(x = Treatment, y = meanAUC, ymin = meanAUC - seAUC, ymax = meanAUC + seAUC,
                   colour = Experiment, group = Experiment))+
  geom_point(size = 5,position = position_dodge(0.25))+
  geom_line(position = position_dodge(0.25))+
  geom_errorbar(width = 0.1, position = position_dodge(0.25))+
  theme_bw(base_size = 25)

# MODEL AUC ----
# Suggestive of Fast - Slow model
aucMod <- lm(AUC ~ Treatment * Experiment, data = AUC)
Anova(aucMod)
heplots::etasq(aucMod, anova = TRUE) #substantial effect sizes

# T- and P-opt values for analysis ----
# use clone_pred

P_T_opts <- pd %>%
  # get the rows where clone preds are max in each group
  filter(clone_pred == max(clone_pred)) %>% 
  select(Temperature, Treatment, Experiment, Clone, clone_pred) %>% 
  # rename to conventions
  rename(Topt = Temperature, Popt = clone_pred)

# Topt/Popt analysis ----
# GENERALIST SPECIALIST Suggestion
mod_Topts <- lm(Topt ~ Treatment * Experiment, data = P_T_opts)
mod_Popts <- lm(Popt ~ Treatment * Experiment, data = P_T_opts)
heplots::etasq(mod_Topts, anova = TRUE) # small effect sizes on Topt
heplots::etasq(mod_Popts, anova = TRUE) # substantial effect sizes Performance

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
  labs(y = expression(paste("Fecundity ", Sigma, "3-clutches")))+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)
