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

# scaleDat <- SGR_Ind %>% as.data.frame() %>% 
#   filter(Clone != "D10_A13") %>% 
#   mutate(maxInduction = scale(maxInduction),
#          Growth0 = scale(Growth0),
#          Growth1 = scale(Growth1),
#          factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
#   na.omit() %>% 
#   data.frame()

scaleDat <- SGR_Ind %>% as.data.frame() %>% 
  filter(Clone != "D10_A13") %>% 
  mutate(factorial = paste(Experiment, Treatment, Clone, sep = ":")) %>% 
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

# expanded temperature range
newX2 <- expand.grid(
  Temperature = seq(5,30,length = 100),
  Treatment = unique(scaleDat$Treatment),
  Experiment = unique(scaleDat$Experiment),
  Clone = unique(scaleDat$Clone)
)

# predicted Fecundities
fixed_pred <- predict(mod, newdata = newX,re.form = NA)
fixed_pred2 <- predict(mod, newdata = newX2, re.form = NA)
clone_pred <- predict(mod, newdata = newX, 
                      re.form = ~(poly(Temperature,2)|Clone))

# housekeeping
pd <- data.frame(newX, fixed_pred, clone_pred) %>% 
  mutate(Experiment = factor(Experiment, levels = c("Acclim", "Acute")))

pd2 <-data.frame(newX2, fixed_pred2, clone_pred) %>% 
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

# align with theory picture
ggplot(pd2, aes(x = Temperature, y = fixed_pred2, colour = Treatment,
                linetype = Experiment))+
  geom_line(size = 2, alpha = 0.3)+
  geom_line(size = 2, data = pd, aes(x = Temperature, y = fixed_pred, colour = Treatment,
                                     linetype = Experiment))+
  geom_vline(xintercept = c(13,28), col = 'grey30')+
  geom_hline(yintercept = 0, col = 'grey30')+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  scale_linetype_manual(values = c(Acute = "dotted", Acclim = "solid"))+
  labs(y =expression(paste("Somatic Growth Rate")))+
  ylim(-0.12,0.12)+
  theme_bw(base_size = 15)

## AUC and Popt/Topt analysis ----

# calcuate area under the curves ----
AUC <- pd %>% group_by(Clone, Treatment, Experiment) %>%
  summarise(AUC = auc(x = Temperature, y = clone_pred))

AUCsum <- AUC %>% group_by(Treatment, Experiment) %>% 
  summarise(meanAUC = mean(AUC),
            seAUC = sd(AUC)/sqrt(sum(!is.na(AUC))))

# GRAPH AUC: predation increases the AUC ----
AUCsum <- AUCsum %>% mutate(Experiment = factor(Experiment, levels = c("Acute","Acclim")))

ggplot(AUCsum, aes(x = Experiment, y = meanAUC, 
                   ymin = meanAUC - seAUC, 
                   ymax = meanAUC + seAUC,
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
summary(aucMod)

# T- and P-opt values for analysis ----
# use clone_pred

# This grabs the peaks for each clone
P_T_opts <- pd %>% group_by(Clone, Treatment, Experiment) %>% 
  # get the rows where clone preds are max in each group
  filter(clone_pred == max(clone_pred)) %>% 
  select(Temperature, Treatment, Experiment, Clone, clone_pred) %>% 
  # rename to conventions
  rename(Topt = Temperature, Popt = clone_pred) %>% 
  select(Treatment, Experiment, Clone, Topt, Popt)

# see the clone specific Popt and Topts
P_T_opts

# This calculates the means and SE for each clone
# of Peak perf and T at which it occurs
P_T_cloneMeans <- P_T_opts %>% 
  group_by(Treatment, Experiment) %>% 
  summarise(
    meanP = mean(Popt),
    seP = sd(Popt)/sqrt(n()),
    meanT = mean(Topt),
    seT = sd(Topt)/sqrt(n())
  )

# Two graphs
Popt_g <- ggplot(P_T_cloneMeans, aes(x = Experiment, y = meanP,
                     ymin = meanP - 1.96*seP,
                     ymax = meanP + 1.96*seP,
                     colour = Treatment, group = Treatment))+
  geom_point(position = position_dodge(width = 0.1))+
  geom_line()+
  geom_errorbar(width = 0, position = position_dodge(width = 0.1))

Topt_g <- ggplot(P_T_cloneMeans, aes(x = Experiment, y = meanT,
                                     ymin = meanT - 1.96*seT,
                                     ymax = meanT + 1.96*seT,
                                     colour = Treatment, group = Treatment))+
  geom_point(position = position_dodge(width = 0.1))+
  geom_line()+
  geom_errorbar(width = 0, position = position_dodge(width = 0.1))

# combined graph
gridExtra::grid.arrange(Popt_g, Topt_g)

# Topt/Popt analysis ----
mod_Topts <- lm(Topt ~ Treatment * Experiment, data = P_T_opts)
mod_Popts <- lm(Popt ~ Treatment * Experiment, data = P_T_opts)
heplots::etasq(mod_Topts, anova = TRUE) # small effect sizes on Topt
heplots::etasq(mod_Popts, anova = TRUE) # small effect sizes Performance

# Theory Plot
ggplot(pd, aes(x = Temperature, y = fixed_pred, 
               group = Treatment, colour = Treatment))+
  # add the curves on top of everything
  geom_line(size = 1)+
  # geom_segment(data = polyg,
  #              aes(x = X1, y = Y1, 
  #              xend = X2, yend = Y1), alpha = 0.3)+
  # geom_segment(data = polyg,
  #              aes(x = X2, y = Y1, 
  #                  xend = X2, yend = Y2), alpha = 0.3)+
      # the rest
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  labs(y = "Somatic Growth Rate (mm/day)")+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)


# Four Traits: Popt, Topt, Min, Max ----

# Gen specialist if Popt up and 13 or 28 up because higher performance at 13 or 28 
# equates with steeper drop and narrower

# get the 13 and 28 values
min_max <- pd %>% filter(Temperature == 13 | Temperature == 28) %>% 
  select(Temperature, Treatment, Experiment, clone_pred) %>% 
  spread(Temperature, clone_pred)

# merge the data frames
four_traits <- left_join(min_max, P_T_opts)

# Test for Gen-Spec: Popt up and 13/28 up; If Popt up and 13/28 same = Fast Slow
# As Popt increases - both 13 and 28 increase... showing G-S
GS1 <- ggplot(four_traits, aes(x = Popt, y = `13`))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Treatment)+
  ggtitle("Gen-Spec 1")

GS2 <- ggplot(four_traits, aes(x = Popt, y = `28`))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Treatment)+
  ggtitle("Gen-Spec 2")

# Test Hot-Cold: Topt up 13/28 down (hits the temp line lower as it's moved over)
# As Topt increases, Performance at 13 decreases... showing fast slow
HC1 <- ggplot(four_traits, aes(x = Topt, y = `13`))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Treatment)+
  ggtitle("Hot-Cold 1")

# As Topt increases, Performance at 28 does not change
HC2 <- ggplot(four_traits, aes(x = Topt, y = `28`))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Experiment ~ Treatment)+
  ggtitle("Hot-Cold 2")

gridExtra::grid.arrange(GS1, GS2, HC1, HC2)

# # T- and P-opt values for plotting (use fixed_pred) ----
# # NOT USING NOW
# graph_PT <- pd %>% ungroup() %>% 
#   select(Temperature, Treatment, Experiment, fixed_pred) %>% 
#   group_by(Treatment, Experiment) %>% 
#   filter(fixed_pred == max(fixed_pred)) %>%
#   rename(Topt = Temperature, Popt = fixed_pred) %>% 
#   distinct()
# 
# # GRAPH Topt/Popt ----
# # a thing of beauty
# 
# ggplot(pd, aes(x = Temperature, y = fixed_pred, 
#                group = Treatment, colour = Treatment))+
#   # add the background clone lines
#   geom_vline(aes(xintercept = Topt, colour = Treatment), data = P_T_opts, 
#              alpha = 0.1)+
#   geom_hline(aes(yintercept = Popt, colour = Treatment), data = P_T_opts,
#              alpha = 0.1)+
#   # add the fixed lines
#   geom_vline(aes(xintercept = Topt, colour = Treatment), data = graph_PT, 
#              alpha = 0.6)+
#   geom_hline(aes(yintercept = Popt, colour = Treatment), data = graph_PT,
#              alpha = 0.6)+
#   # add the curves on top of everything
#   geom_line(size = 1)+
#   # the rest
#   scale_colour_manual(values = c(Control = "black", Predator = "red"))+
#   labs(y = "Somatic Growth Rate (mm/day)")+
#   facet_wrap(~Experiment, ncol = 2)+
#   theme_bw(base_size = 15)
# 
# # polyg <- pd %>% group_by(Treatment, Experiment) %>% 
# #   summarise(
# #     X1 = min(Temperature), X2 = max(Temperature),
# #     Y1 = min(fixed_pred), Y2 = last(fixed_pred)
# #   )
