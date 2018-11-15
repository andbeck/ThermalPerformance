# Hebe All Data (18.7.18)
#rm(list=ls())
library(lme4)
library(splines)
library(tidyverse)
library(car)


# all data ---------------------------------------------------------------
allDat <- read_csv("AllResults_18_07_18.csv")
#allDat <- read_csv(file.choose())

# check your data and levels of the factors ---------------------------------------------------------------
# distinct(allDat, Clone)
# distinct(allDat, Treatment)
# distinct(allDat, Temperature)
# distinct(allDat, Experiment)
# distinct(allDat, Stage) # problems
# distinct(allDat, Mature) # problems

# fix the factor levels and get rid of D10_74 ---------------------------------------------------------------
allDat2 <- 
  allDat %>%
  mutate(Experiment = ifelse(Experiment == "Aclim", "Acclim", Experiment)) %>% 
  mutate(Stage = ifelse(Stage == "Mature", "mt", Stage)) %>% 
  mutate(Clone = ifelse(Clone == "LD34"| Clone == "LD35", "LD33", Clone)) %>% 
  filter(Clone != "D10_74")

# Isolate Maturation Data (Size/Age) ---------------------------------------------------------------
matDat <- filter(allDat2, Stage == "mt") %>% as.data.frame()
matDat

#### VERSION 1: Clone Means Approach ---------------------------------------------------------------
# This is a simpler way.  We assume each CLONE has a mean
# this gives ~ 10 data points per temp*treat*exp combination
cloneMeans <- matDat %>% 
  group_by(Clone, Treatment, Temperature, Experiment) %>% 
  summarise(
    meanSize = mean(Body, na.rm = TRUE),
    meanAge = mean(Age, na.rm = TRUE))

# these are simpler to understand than clone specific
SaM_clone <- ggplot(cloneMeans, aes(x = Temperature, y = meanSize, colour = Treatment))+
  geom_point()+
  geom_smooth(method = gam, formula = y ~ s(x, k = 3)) +
  #geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE)+
  facet_wrap(~Experiment)

AaM_clone <- ggplot(cloneMeans, aes(x = Temperature, y = meanAge, colour = Treatment))+
  geom_point()+
  geom_smooth(method = gam, formula = y ~ s(x, k = 3)) +
  #geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE)+
  facet_wrap(~Experiment)

gridExtra::grid.arrange(SaM_clone, AaM_clone)

# Model Clone Means ---------------------------------------------------------------
# does not require mixed effects
# fitting polynomial models to the four panel graph.

modSize_clone <- lm(meanSize ~ poly(Temperature,2)* Treatment*Experiment, 
                    data = cloneMeans)

# Temeprature effect is non-linear
# Treatment changes mean Size
# Experiment changes mean Size
# all effects additive
Anova(modSize_clone)

modAge_clone <- lm(meanAge ~ poly(Temperature,2)* Treatment*Experiment, data = cloneMeans)

# Temperature effect is non-linear
# marginal change in mean Age by experiment
# NO effect of Treatment
Anova(modAge_clone)

#### VERSION 2: Clonal Variation approach (mixed effects) ---------------------------------------------------------------
# Raw Data and smoothing fits (polynomial) 
# paints a rather complex picture, especially for Size
# omitting geom_point() just makes a cool and clean picture!
SaM <- ggplot(matDat, aes(x = Temperature, y = Body, 
                   colour = Treatment))+
  #geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()

AaM <- ggplot(matDat, aes(x = Temperature, y = Age, 
                          colour = Clone))+
  #geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Treatment~Experiment)+
  theme_bw()

gridExtra::grid.arrange(SaM, AaM)

# Toy Models ---------------------------------------------------------------
# temperature as a polynomial interacting with treatment and experiment
# each clone has it's own average response to temperature
# This estiamtes genetic variation and is a complex model

# Age Model
modAge_poly <- lmer(Age ~ poly(Temperature, 2)*Treatment*Experiment + 
                      (poly(Temperature, 2)|Clone), 
                    data = matDat)

## polynomial justified
#modAge_lin <- lmer(Age ~ Temperature*Treatment*Experiment + (poly(Temperature, 2)|Clone), data = matDat)
#anova(modAge_poly, modAge_lin)

# Size Model
modSize_poly <- lmer(Body ~ poly(Temperature, 2)*Treatment*Experiment + 
                       (poly(Temperature, 2)|Clone), 
                     data = matDat)




# fits <- matDat %>% group_by(Clone) %>% 
#   do(fitTemps = lm(Body ~ Temperature, data = .))

## polynomial justified
# modSize_lin <- lmer(Body ~ Temperature*Treatment*Experiment + (poly(Temperature, 2)|Clone), data = matDat)
# anova(modSize_poly, modSize_lin)

# Age Effects (very similar to Version 1)
# ** effect of Temp is nonlinear and varies by Treatment.  
# Experiments produce different average Age
Anova(modAge_poly, test.statistic = "F")

# Size Efects (New Treatment x Temeprature effect here compared to Version 1)
# ** Temperature effect is non-linear and varies by Treatment BUT NOT experiment
# ** mean size varies between Experiment
Anova(modSize_poly, test.statistic = "F")


# spline mixed approach
mod_spline <- lmer(Body ~ bs(Temperature)*Treatment + 
                     Experiment +
                     (bs(Temperature)|Clone),
                   data = matDat)

summary(mod_spline)
Anova(mod_spline)
ranef(mod_spline)
plot(mod_spline)

newX <- expand.grid(
  Temperature = seq(13,28,1),
  Treatment = unique(matDat$Treatment),
  Experiment = unique(matDat$Experiment),
  Clone = unique(matDat$Clone)
)

fixed_pred <- predict(mod_spline, newdata = newX, re.form = NA)
clone_pred <- predict(mod_spline, newdata = newX, re.form = ~(bs(Temperature)|Clone))

pd <- data.frame(newX, fixed_pred, clone_pred)

ggplot(pd, aes(x = Temperature, y = fixed_pred))+
  geom_line(size = 2)+
  geom_line(aes(x = Temperature, y = clone_pred, colour = Clone), 
            size = 1, alpha = 0.5)+
  facet_grid(Experiment ~ Treatment)

# AUCs
AUC <- pd %>% group_by(Treatment, Experiment) %>% 
  summarise(AUC = auc(x = Temperature, y = clone_pred))

ggplot(AUC, aes(x = Treatment, y = AUC, 
                colour = Clone, group = Clone))+
  geom_point()+
  geom_line()+
  facet_wrap(~ Experiment)

