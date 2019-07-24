library(tidyverse)
library(lme4)
library(car)
library(heplots)
library(MESS)
library(lmerTest)

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")
names(SGR_Ind)

# Create a scaled dataset and a response one with Factorial ----------------------------
# Omit D10_A13 because so much missing

scaleInd_Pred <- SGR_Ind %>% 
  arrange(Clone, Temperature, Treatment, Experiment) %>% 
  mutate(ID = paste(Clone, Temperature, Treatment, Experiment, sep = "-")) %>% 
  select(ID, Clone, Temperature, Treatment, Experiment, maxInduction) %>% 
  filter(Treatment == "Predator")

# Fec Data
source("./Scripts/Working/MakeFecData_APB.R")

# consider using Fec or Effort... both give same result - no trade-off.

Fec_scale_pred <- Fec_scale %>% ungroup() %>% 
  select(Clone, Temperature, Treatment, Experiment, Fec, Effort) %>% 
  arrange(Clone, Temperature, Treatment, Experiment) %>% 
  mutate(ID = paste(Clone, Temperature, Treatment, Experiment, sep = "-")) %>% 
  filter(Treatment == "Predator")

glimpse(Fec_scale_pred)

ggplot(Fec_scale_pred, aes(x = Temperature, y = Effort))+
  geom_jitter(size = 2, width = 0.5, aes(colour = Clone), alpha = 0.3)+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, colour = "black")+
  geom_line(stat = 'smooth', aes(group = Clone, colour = Clone),
            method = lm, formula = y ~ poly(x, 2), se = FALSE, alpha = 0.3)+
  facet_grid(~Experiment)+
  theme_bw(base_size = 10)+
  guides(colour = FALSE)


###
# look at the data
glimpse(scaleInd_Pred)
glimpse(Fec_scale_pred)

tempInd <- filter(scaleInd_Pred, 
       Clone == "LD33", 
       Experiment == "Acute", 
       Temperature == 13) %>% data.frame()

tempFec <- filter(Fec_scale_pred, 
       Clone == "LD33", 
       Experiment == "Acute", 
       Temperature == 13) %>% data.frame()

tempFec
tempInd

# CHECK THIS AGAIN.
fec_ind <- full_join(scaleInd_Pred, Fec_scale_pred)
glimpse(fec_ind)

# test
filter(fec_ind, 
       Clone == "LD33", Treatment == "Predator", 
       Experiment == "Acute", Temperature == 13)

# plot all of the data
ggplot(fec_ind, aes(x = maxInduction, y = Effort))+
  geom_smooth(method = lm, se = FALSE)+
  geom_jitter(aes(colour = Clone))+
  facet_grid(Experiment~Temperature)

# take clone means
mean_fec_ind <- fec_ind %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  summarise(meanInd = mean(maxInduction, na.rm = TRUE),
            meanFec = mean(Effort, na.rm = TRUE))

# factorise by temperature
F0 <-  mean_fec_ind %>% 
  ungroup() %>% 
  mutate(Temperature = factor(Temperature)) %>% 
  filter(Treatment == "Predator")

# plot it
# perhaps a strong ind x temp in acclim, but not in acute?
ggplot(F0, aes(x = meanInd, y = meanFec))+
  # coloured points by temperature
  geom_point(aes(group = Temperature, colour = Temperature))+
  # temperature regressions among clones
  geom_smooth(aes(group = Temperature, colour = Temperature), 
              method = lm, se=FALSE, alpha = 0.5)+
  # overall regression
  geom_smooth(method = lm, se = FALSE, col = 'black')+
  # customise colour scheme
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(Experiment ~ .)

# Models of clone means
modTO <- lm(meanFec ~ meanInd * Temperature * Experiment, data = F0)

Anova(modTO)
summary(modTO)

newX <- expand.grid(
  meanInd = seq(from = 0, to = 100, by = 10),
  Temperature = unique(F0$Temperature),
  Experiment =  unique(F0$Experiment))

newY <- predict(modTO, newdata = newX, interval = 'confidence')

plotThese <- data.frame(newX, newY) %>% 
  rename(predictedFec = fit, Induction = meanInd)

ggplot(plotThese, aes(x = Induction, y = predictedFec, 
                      group = Temperature, colour = Temperature))+
  geom_line()+
  geom_point(data = F0, 
             aes(x = meanInd, y = meanFec))+
  facet_wrap(~Experiment)

# lmer model ---------------------------------------------------------
modTradeOff <- lmer(Fec ~ maxInduction*Temperature*Experiment+
                      (maxInduction|Clone), data = fec_ind, na.action = 'na.omit',
                    control = lmerControl(optimizer = "Nelder_Mead"))
modTradeOff2 <- lmer(Fec ~ (maxInduction+Temperature+Experiment)^2+
                      (maxInduction|Clone), data = fec_ind,na.action = 'na.omit',
                    control = lmerControl(optimizer = "Nelder_Mead"))

Anova(modTradeOff, test ="F")
anova(modTradeOff, modTradeOff2, test = "Chisq")

# use lmerTest model
mod3 <- lmerTest::lmer(Fec ~ maxInduction*Temperature*Experiment+
                         (maxInduction|Clone), data = fec_ind, na.action = 'na.omit',
                       control = lmerControl(optimizer = "Nelder_Mead"))

anova(mod3, type = "II")


# prep plot lmer  Result 

newX <- expand.grid(
  maxInduction = seq(from = 0,
                     to = 100,
                     length = 10),
  Temperature = unique(fec_ind$Temperature),
  Experiment = unique(fec_ind$Experiment),
  Clone = unique(fec_ind$Clone))

fixed_pred <- predict(modTradeOff, newdata = newX, re.form = NA)
clone_pred <- predict(modTradeOff, newdata = newX, 
                      re.form = ~(maxInduction|Clone))

pd <- data.frame(newX, fixed_pred, clone_pred) %>% 
  mutate(Experiment = factor(Experiment, levels = c("Acclim","Acute")))

# lmer model Plot 
lmerFixed <- ggplot(pd, aes(x = maxInduction, y = fixed_pred, 
                            colour = factor(Temperature), 
                            group = Temperature))+
  geom_line(size = 2)+
  geom_point(data = fec_ind, aes(x = maxInduction, y = Fec), 
             colour = "grey")+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  #scale_y_continuous(breaks = seq(from = -2, to = 2, by = 0.5))+
  facet_grid(Experiment ~ Temperature)+
  labs(y = "Reproductive Effort", x = "Max Induction")+
  theme_bw(base_size = 15)+
  guides(color=guide_legend(title="Temp C˚"))+
  theme(legend.position = "top")


lmerFixed

lmerClone <- ggplot(pd, aes(x = maxInduction, y = clone_pred, 
                            colour = factor(Temperature), 
                            group = Temperature))+
  geom_line(size = 2)+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(factor(Experiment, levels = c("Acclim","Acute")) ~ Clone)+
  labs(y = "Reproductive Effort", x = "Max Induction")+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

lmerClone

gridExtra::grid.arrange(lmerFixed, lmerClone, ncol = 1)

# MCMCglmm model --------------------------------
# 4 terms in the random effect (intercept + 3)
prior <- list(R = list(V = 1, n = 0.002),
              G = list(G1 = list(V = diag(1), n = 1)))

modBay <- MCMCglmm(Fec ~ maxInduction*Temperature*Experiment,
                   random = ~us(maxInduction):Clone, 
                   family = 'gaussian', prior = prior, pr = TRUE,
                   data = data.frame(na.omit(fec_ind)))

modBay2 <- MCMCglmm(Fec ~ (maxInduction+Temperature+Experiment)^2,
                    random = ~us(maxInduction):Clone, 
                    family = 'gaussian', prior = prior, pr = TRUE,
                    data = data.frame(na.omit(fec_ind)))

summary(modBay)$sol
modBay$DIC-modBay2$DIC

## Fecundity Plasticity vs. Morphology ----
# trying to get the change in effort plotted against the induced morphology
# do clones that are more plastic in reproduction have lower morphology

AccCon <- Fec_scale %>% 
  filter(Treatment == "Control", Experiment == "Acclim") %>% 
  select(Effort) %>% 
  ungroup() %>% 
  select(-ID)

AccPred <- Fec_scale %>% 
  filter(Treatment == "Predator", Experiment == "Acclim") %>% 
  select(Effort) %>% 
  rename(EffortP = Effort) %>% 
  ungroup() %>% 
  select(-ID)

AccPredMorph <- scaleInd_Pred %>% 
  filter(Experiment == "Acclim")

AccCon$ID
AccPred$ID

out <- full_join(AccCon, AccPred, by = c("Clone","Temperature")) %>% 
  mutate(deltaEffort = Effort - EffortP)

out2 <- full_join(out, AccPredMorph, by = c("Clone", "Temperature"))

outPlot <- out2 %>% 
  group_by(Clone, Temperature) %>% 
  summarise(
    meanDE = mean(deltaEffort, na.rm = TRUE),
    meanInd = mean(maxInduction, na.rm = TRUE)
  )

# WHAT do we expect the plasticity to be?  ±?
ggplot(outPlot, aes(x = meanInd, y = meanDE))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  theme_bw(base_size = 15)

####

##Build Stearns like picture with 13 vs. 28 C ----
stearns <- fec_ind %>% filter(Temperature == 13|Temperature == 28)
glimpse(stearns)

ggplot(stearns, aes(x = maxInduction, y = Effort, colour = factor(Temperature)))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE)+
  facet_grid(~Experiment)

stearns_predictions <- pd %>% filter(Temperature == 13|Temperature == 28)
glimpse(stearns_predictions)

ggplot(stearns_predictions, aes(x = maxInduction, y = fixed_pred, colour = Experiment))+
  geom_smooth(method = 'lm', se = FALSE)+
  facet_grid(Experiment~Temperature)
