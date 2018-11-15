# hebe 

library(tidyverse)
library(splines)
library(lme4)
library(car)
library(MESS)

# data
Ro <- read_csv("./Data/Ro_Data.csv")
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")

## Pop Growth Rate - Ro -----------------------------------------------------
# graphs to start
Ro_3 <- ggplot(Ro, aes(x = Temperature, y = Ro.B3, colour = Treatment,
                       group = Treatment))+
  geom_point() + ylim(0,50) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x,2))+
  facet_wrap( ~ Experiment)+
  theme_bw()


Ro_2 <- ggplot(Ro, aes(x = Temperature, y = Ro.B2, colour = Treatment,
                       group = Treatment))+
  geom_point() + ylim(0,50) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x,2))+
  facet_grid(Experiment ~ Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  ylab("Population Growth Rate")+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45))
  
Ro_2

Ro_1 <- ggplot(Ro, aes(x = Temperature, y = Ro.B1, colour = Treatment,
                       group = Treatment))+
  geom_point() + ylim(0,50) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x,2))+
  facet_wrap( ~ Experiment)+
  theme_bw()


gridExtra::grid.arrange(Ro_1, Ro_2, Ro_3, ncol = 3)


# model Ro ------------------------------------------------------------

# CHOOSE ONE
# model with spline
mod_Ro2 <- lmer(Ro.B2 ~ poly(Temperature, 2) * Treatment * Experiment +
                  (poly(Temperature,2)+Treatment+Experiment|Clone),
              data = Ro)

par(mfrow = c(2,2))
plot(mod_Ro2)

Anova(mod_Ro2, test.statistic = "F")

# plot Ro Model --------------------

newX <- expand.grid(
  Temperature = seq(13,28,length = 100),
  Treatment = unique(Ro$Treatment),
  Experiment = unique(Ro$Experiment),
  Clone = unique(Ro$Clone)
)

# predicted Ro ----
fixed_pred <- predict(mod_Ro2, newdata = newX,re.form = NA)
clone_pred <- predict(mod_Ro2, newdata = newX, re.form = ~(Experiment+Treatment+poly(Temperature,2)|Clone))


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



# nm <- function(Temperature, p, Tmax, dT, LL){ exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT)+LL}
# Temp <- seq(14, 28, 1)
# df2 <- filter(Ro, Experiment == "Acclim"&Treatment == "Control")
# 
# df <- data.frame(Temp, Perf = nm(Temp, p = 0.01, Tmax = 25, dT = 2, LL = 10))
# ggplot(df, aes(x = Temp, y = Perf))+
#   geom_line()+
#   geom_point(data = df2, aes(x = Temperature, y = Ro.B2))
#   
# 
# nls_Ro2 <- nls(Ro.B2~ exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT)+LL, 
#                start = list(p = 0.14, Tmax = 28, dT = 4, LL = 5), 
#                data = filter(Ro, Experiment == "Acclim"&Treatment == "Control"),
#                trace = TRUE)


# library(devRate)
# 
# myTemp = select(df2, Temperature) %>% as.numeric()
# myDev = select(df2, Ro.B2) %>% as.numeric()
# devRateModel(eq = campbell_74, temp = myT, devRate = myDev)


# model via polynomial
# mod_Ro2 <- lm(Ro.B2 ~ poly(Temperature,2) * Treatment * Experiment,
#               data = Ro)


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

ggplot(pd, aes(x = Temperature, y = fixed_pred, group = Treatment))+
  geom_line(size = 2)+
  geom_vline(aes(xintercept = Temperature), data = Topts, col = 'red')+
  geom_hline(aes(yintercept = fixed_pred), data = Topts, col = 'blue')+
  labs(y = "Predicted Population Growth Rate (mm/day)")+
  facet_wrap(~Experiment, ncol = 2)+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")



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

# second show how the optimal temperature (where the peak happens in temperature space)
# changes among the treatments...
# NOTE how the optimal temperature for Ro is at the min Temp for both predator and no predator under
# acute, but that it's higher under acclim, and that predation lowers the temperature at which max Ro 
# occurs!

TPC2 <- ggplot(filter(out, Metric == "Opt"), 
       aes(x = Treatment, y = Temperature, colour = Experiment, group = Experiment))+
  geom_point()+geom_line()+
  ylab("Temperature at which Max Ro Occurs (Topt)")+
  theme_bw()

gridExtra::grid.arrange(TPC1, TPC2, ncol = 2)

# ALL RESULTS 1 PANE
gridExtra::grid.arrange(Results_Plot, AUC_plot, TPC1, TPC2, ncol = 2)


## SGR ------------------------------------------

SGR_finSize <- ggplot(SGR_Ind, aes(x = Temperature, y = Growth0, 
                                   colour = Clone))+
  geom_point()+
  geom_smooth(span = 1, se = FALSE)+
  facet_grid(Treatment ~ Experiment+Clone)+
  theme_bw()

SGR_instar1 <- ggplot(SGR_Ind, aes(x = Temperature, y = Growth1, 
                           colour = Clone))+
  geom_point()+
  geom_smooth(span = 1, se = FALSE)+
  facet_grid(Treatment ~ Experiment+Clone)+
  theme_bw()

SGR_instar2 <- ggplot(SGR_Ind, aes(x = Temperature, y = Growth2, 
                                   colour = Clone))+
  geom_point()+
  geom_smooth(span = 1, se = FALSE)+
  facet_grid(Treatment ~ Experiment+Clone)+
  theme_bw()

gridExtra::grid.arrange(SGR_finSize, SGR_instar1, SGR_instar2)


ggplot(SGR_Ind, aes(x = Temperature, y = Growth1, colour = Treatment,
                    group = Treatment))+
  geom_point()+
  geom_smooth(method = lm, se=F)+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(14,16,20,24,28))+
  ylab("Somatic Growth Rate")+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45))


# Model SGR
model_SGR <- lmer(Growth2 ~ bs(Temperature)*Treatment*Experiment+
                    (Experiment + Treatment + bs(Temperature)|Clone), data = SGR_Ind)


# Anova from car and the F-test via kenward-rogers (everyone expects this)
Anova(model_SGR, test.statistic = "F")

summary(model_SGR)


# fit model with simplest random effect justified by designe
# (1|Clone)
model_SGR_Simple <- lmer(Growth2 ~ bs(Temperature)*Treatment*Experiment+
                         (1|Clone), data = SGR_Ind)

# compare complex to simple
anova(model_SGR, model_SGR_Simple)

# observed vs. predicted
plot(predict(model_SGR, type = 'response') ~ na.omit(SGR_Ind$Growth2))
abline(0,1)

# plot SGR
newX <- expand.grid(
  Temperature = seq(13,28,length = 50),
  Treatment = unique(SGR_Ind$Treatment),
  Experiment = unique(SGR_Ind$Experiment),
  Clone = unique(SGR_Ind$Clone)
)

# new Y's - one for average (fixed) and one for each clone (re.form = ~....)
fixed_pred <- predict(model_SGR, newdata = newX, re.form = NA)
clone_pred <- predict(model_SGR, newdata = newX, re.form = ~(Experiment+Treatment+bs(Temperature)|Clone))

# housekeeping
pd <- data.frame(newX, fixed_pred, clone_pred)

ggplot(pd, aes(x = Temperature, y = fixed_pred))+
  geom_line(size = 2)+
  #geom_line(aes(x = Temperature, y = clone_pred, colour = Clone), 
  #          size = 1, alpha = 0.5)+
  labs(y = "Predicted Somatic Growth Rate (mm/day)")+
  facet_grid(Experiment ~ Treatment)+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")



# AUCs
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

ggplot(pd, aes(x = Temperature, y = fixed_pred, group = Treatment))+
  geom_line(size = 2)+
  geom_vline(aes(xintercept = Temperature), data = Topts, col = 'red')+
  geom_hline(aes(yintercept = fixed_pred), data = Topts, col = 'blue')+
  labs(y = "Predicted Somatic Growth Rate (mm/day)")+
  facet_wrap(~Experiment, ncol = 4)+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")



# Induction ------------------------------------------

ggplot(SGR_Ind, aes(x = Temperature, y = maxInduction, 
                           colour = Clone))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(Treatment ~ Experiment+Clone)+
  theme_bw()

# Model Ind
model_Ind <- lmer(maxInduction ~ Temperature*Treatment*Experiment+
                    (Experiment+Treatment+Temperature|Clone), data = SGR_Ind)

Anova(model_Ind, test.statistic = "F")

summary(model_Ind)

# test random effects
model_Ind_Simple <- lmer(maxInduction ~ Temperature*Treatment*Experiment+
                         (1|Clone), data = SGR_Ind)

anova(model_Ind, model_Ind_Simple)

# plot IND
newX <- expand.grid(
  Temperature = seq(13,28,1),
  Treatment = unique(SGR_Ind$Treatment),
  Experiment = unique(SGR_Ind$Experiment),
  Clone = unique(SGR_Ind$Clone)
)

fixed_pred <- predict(model_Ind, newdata = newX, re.form = NA)
clone_pred <- predict(model_Ind, newdata = newX, re.form = ~(Experiment+Treatment+Temperature|Clone))

pd <- data.frame(newX, fixed_pred, clone_pred)

ggplot(pd, aes(x = Temperature, y = fixed_pred))+
  geom_line(size = 2)+
  geom_line(aes(x = Temperature, y = clone_pred, colour = Clone), 
            size = 1, alpha = 0.5)+
  facet_grid(Treatment~Experiment)+
  labs(y = "Maximum Induction (Defence Morphology) Score")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")


# Trade-off
pred_only <- filter(SGR_Ind, Treatment == "Predator")

# lots of variation in relationship
# including reversals
ggplot(pred_only, aes(x = maxInduction, y = Growth2, colour = Clone))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_grid(Experiment ~ Clone)

# what's up with temperature
ggplot(pred_only, aes(x = maxInduction, y= Growth2))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_grid(Experiment ~ Temperature)

# more what's up with temperature
ggplot(pred_only, aes(x = maxInduction, y= Growth2, colour = Temperature,
                      group = Temperature))+
  geom_point(size = 5)+
  geom_smooth(method = lm, se = FALSE)+
  facet_grid(~Experiment)


# model
modTrad <- lmer(Growth2 ~ maxInduction*Experiment +
                  (Experiment+maxInduction|Clone), data = SGR_Ind)


# temperature effects
modTrad2 <- lmer(Growth2 ~ maxInduction*Temperature*Experiment+
                   (maxInduction+Temperature|Clone), data = pred_only)

Anova(modTrad2)
summary(modTrad2)

# tests
Anova(modTrad, test.statistic = "F")
summary(modTrad)

# plot TradeOff Result
newX <- expand.grid(
  maxInduction = seq(0,100, 1),
  Experiment = unique(SGR_Ind$Experiment),
  Temperature = seq(14,28,2),
  Clone = unique(SGR_Ind$Clone)
)

fixed_pred <- predict(modTrad, newdata = newX, re.form = NA)
clone_pred <- predict(modTrad, newdata = newX, re.form = ~(Experiment+maxInduction|Clone))

# with temperature effects
fixed_pred <- predict(modTrad2, newdata = newX, re.form = NA)
clone_pred <- predict(modTrad2, newdata = newX, re.form = ~(1|Clone))

pd <- data.frame(newX, fixed_pred, clone_pred)

ggplot(pd, aes(x = maxInduction, y = fixed_pred))+
  geom_line(size = 2)+
  geom_line(aes(x = maxInduction, y = clone_pred, colour = Clone), 
            size = 1, alpha = 0.5)+
  facet_grid(~Experiment)+
  labs(y = "Ro")+
  theme_bw()


ggplot(pd, aes(x = maxInduction, y = fixed_pred, colour = Temperature))+
  geom_point()+
  facet_wrap(~Experiment)+
  theme_bw()
