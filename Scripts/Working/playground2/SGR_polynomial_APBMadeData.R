# Somatic Growth Rate
# log(size difference)/age difference

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

## SGR Analyses ----

# look at the data
glimpse(SGR_scale)
p_Gr # the plot

# The model ----

# Effect of Temp varies by treatment and by Experiment
# Effect of Treatment does not vary by Experiment (Predation effect is invariant)
mod <- lmer(Growth0 ~ (poly(Temperature,2)+Treatment+Experiment)^3 + 
              (poly(Temperature,2)|Clone), 
            data = SGR_scale)

summary(mod)  
Anova(mod)

# Plot the results ----

# showing variation among clones
# and fixed effect

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(SGR_scale$Treatment),
  Experiment = unique(SGR_scale$Experiment),
  Clone = unique(SGR_scale$Clone)
)

# predicted Fecundities
fixed_pred <- predict(mod, newdata = newX, re.form = NA)
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
  facet_grid(Experiment ~ Treatment)+
  labs(y ="Somatic Growth Rate (mm/day)")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")

