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
library(heplots)

## Source the Trait making script ----

# this will make the fecundity data
# it will produce a plot too
# takes a few seconds.

source("./Scripts/Working/MakeFecData_APB.R")

## Fecundity Analyses ----

# look at the data
glimpse(Fec_scale)
p_Fec # the plot

# USING Fec OR Fec/Body Size for Reproductive Effort ----
# Reproductive Effort gives substantial three way result

ggplot(Fec_scale, aes(x = Temperature, y = Fec/Body, 
                     colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  ylab(expression(paste("Fecundity (", Sigma, "3-clutches)")))+
  theme_bw()

# The model ----

# Effect of Temp varies by treatment and by Experiment
# Effect of Treatment does not vary by Experiment (Predation effect is invariant)
mod <- lmer(Fec/Body ~ (poly(Temperature,2)+Treatment+Experiment)^3 + 
              (poly(Temperature,2)|Clone), 
            data = Fec_scale, control = lmerControl(optimizer = "Nelder_Mead"))

summary(mod)  
Anova(mod, test = "F")

# Plot the results ----

# showing variation among clones
# and fixed effect

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(Fec_scale$Treatment),
  Experiment = unique(Fec_scale$Experiment),
  Clone = unique(Fec_scale$Clone)
)

# use with theory picture?
# expanded temperature range
newX2 <- expand.grid(
  Temperature = seq(5,30,length = 100),
  Treatment = unique(Fec_scale$Treatment),
  Experiment = unique(Fec_scale$Experiment),
  Clone = unique(Fec_scale$Clone)
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
  labs(y =expression(paste("Fecundity (", Sigma, "3-clutches)")))+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

# align with theory picture (uses wider range of the Temperatures)
ggplot(pd2, aes(x = Temperature, y = fixed_pred2, colour = Treatment,
               linetype = Experiment))+
  geom_line(size = 2, alpha = 0.3)+
  geom_line(size = 2, data = pd, aes(x = Temperature, y = fixed_pred, colour = Treatment,
                                     linetype = Experiment))+
  geom_vline(xintercept = c(13,28), col = 'grey30')+
  geom_vline(xintercept = c(15,24), col = 'grey30')+
  geom_hline(yintercept = c(0,8), col = 'grey30')+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  scale_linetype_manual(values = c(Acute = "dashed", Acclim = "solid"))+
  labs(y =expression(paste("Fecundity (", Sigma, "3-clutches)")))+
  ylim(-1,30)+
  theme_bw(base_size = 15)


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
  summarise(AUC = auc(x = Temperature, y = clone_pred))

AUCsum <- AUC %>% group_by(Treatment, Experiment) %>% 
  summarise(meanAUC = mean(AUC),
            seAUC = sd(AUC)/sqrt(sum(!is.na(AUC))))

# GRAPH AUC: predation increases the AUC ----
AUCsum <- AUCsum %>% mutate(Experiment = factor(Experiment, levels = c("Acute","Acclim")))

ggplot(AUCsum, aes(x = Experiment, y = meanAUC, ymin = meanAUC - seAUC, ymax = meanAUC + seAUC,
                   colour = Treatment, group = Treatment))+
  geom_point(size = 5,position = position_dodge(0.25))+
  geom_line(position = position_dodge(0.25))+
  geom_errorbar(width = 0.1, position = position_dodge(0.25))+
  scale_colour_manual(values = c(Control = "black", Predator = "red"))+
  ylab("AUC Fecundity")+
  theme_bw(base_size = 15)

# MODEL AUC ----
# Suggestive of Fast - Slow model
aucMod <- lm(AUC ~ Treatment * Experiment, data = AUC)
Anova(aucMod)
heplots::etasq(aucMod, anova = TRUE) #substantial effect sizes

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
P_T_across_cloneMeans <- P_T_opts %>% 
  group_by(Treatment, Experiment) %>% 
  summarise(
    meanP = mean(Popt),
    seP = sd(Popt)/sqrt(n()),
    meanT = mean(Topt),
    seT = sd(Topt)/sqrt(n())
  )

# Topt/Popt analysis ----
# GENERALIST SPECIALIST Suggestion
mod_Topts <- lm(Topt ~ Treatment * Experiment, data = P_T_opts)
mod_Popts <- lm(Popt ~ Treatment * Experiment, data = P_T_opts)
heplots::etasq(mod_Topts, anova = TRUE) # small effect sizes on Topt
heplots::etasq(mod_Popts, anova = TRUE) # substantial effect sizes Performance

# Four Traits: Popt, Topt, Min, Max ----

# Gen specialist if Popt up and 13 or 28 up because higher performance at 13 or 28 
# equates with steeper drop and narrower

# get the 13 and 28 values
min_max <- pd %>% filter(Temperature == 13| Temperature == 28) %>% 
  select(Temperature, Treatment, Experiment, clone_pred) %>% 
  spread(Temperature, clone_pred)

# merge the data frames
four_traits <- left_join(min_max, P_T_opts)

# get p13/p28 and Popt/Topt
template <- four_traits %>% 
  group_by(Treatment, Experiment) %>% 
  summarise(
    p13 = mean(`13`),
    p28 = mean(`28`),
    meanT = mean(Topt),
    meanP = mean(Popt)
  )

# GS/FS/HC testing
g1 <- ggplot(template, aes(x = meanP, y = p13, colour = Treatment, 
                       shape = Experiment, group = Experiment))+
    geom_point(size = 5)+
    geom_line(colour = 'black')

g2 <- ggplot(template, aes(x = meanP, y = p28, colour = Treatment, 
                           shape = Experiment, group = Experiment))+
  geom_point(size = 5)+
  geom_line(colour = 'black')

gridExtra::grid.arrange(g1,g2)


# # Theory Plot
# ggplot(pd, aes(x = Temperature, y = fixed_pred, 
#                group = Treatment, colour = Treatment))+
#   # add the curves on top of everything
#   geom_line(size = 1)+
#   # geom_segment(data = polyg,
#   #              aes(x = X1, y = Y1, 
#   #                  xend = X2, yend = Y1), alpha = 0.3)+
#   # geom_segment(data = polyg,
#   #              aes(x = X2, y = Y1, 
#   #                  xend = X2, yend = Y2), alpha = 0.3)+
#   # geom_segment(data = polyg,
#   #              aes(x = X1, y = Y1, 
#   #                  xend = X1, yend = Y3), alpha = 0.3)+
#   # the rest
#   scale_colour_manual(values = c(Control = "black", Predator = "red"))+
#   labs(y = expression(paste("Fecundity ", Sigma, "2-clutches")))+
#   facet_wrap(~Experiment, ncol = 2)+
#   theme_bw(base_size = 15)
# 
# 
# # T- and P-opt values for plotting (use fixed_pred) ----
# # NOT USING
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
#     # add the background clone lines
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
#   labs(y = expression(paste("Fecundity ", Sigma, "3-clutches")))+
#   facet_wrap(~Experiment, ncol = 2)+
#   theme_bw(base_size = 15)
# 
# # polyg <- pd %>% group_by(Treatment, Experiment) %>% 
# #   summarise(
# #     X1 = min(Temperature), X2 = max(Temperature),
# #     Y1 = min(fixed_pred), Y2 = last(fixed_pred),
# #     Y3 = first(fixed_pred)
# #   )
