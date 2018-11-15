# Project data analyisis 07/07/18

# libraries we always use
library(tidyverse)
library(ggfortify) 
library(lme4)
library(car)

data <- read.csv("AllResults_16_07_18.csv")
data1 <- data %>% filter(!is.na(Stage))

data1$Stage

maturityData <- data1 %>% filter(Stage == "mt"|Stage == "Mature")
maturityData$Body = as.numeric(as.character(maturityData$Body))
maturityData$Age = as.numeric(as.character(maturityData$Age))
maturityData$Temperature <- factor(maturityData$Temperature)

head(maturityData)

# raw data
ggplot(maturityData, aes(x= Temperature, y=Body,
                             group = Treatment, colour = Treatment)) +
  geom_point() +
  facet_grid(Experiment ~ Clone) +
  ylab("Size at Maturity (mm)")

# size at maturity
maturityDataMean <- maturityData %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanBody = mean(Body, na.rm = TRUE),
            seBody = sd(Body, na.rm = TRUE)/sqrt(sum(!is.na(Body))))
maturityDataMean

ggplot(maturityDataMean, aes(x= Temperature, y= meanBody,
                             group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
       geom_errorbar(aes(ymin = meanBody - seBody, ymax = meanBody + seBody), 
                     width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Size at Maturity (mm)")

# model sizemat
sizeAccute <- filter(maturityData, Experiment == 'Acute')
sizeAcclim <- filter(maturityData, Experiment == "Acclim")

sizeMod_accute <- lmer(Body ~ Temperature * Treatment + (1|Clone), data = sizeAccute)
sizeMod_acclim <- lmer(Body ~ Temperature * Treatment + (1|Clone), data = sizeAcclim)
Anova(sizeMod_accute, test = "F")
Anova(sizeMod_acclim, test = "F")

### age at maturity ###
maturityAgeMean <- maturityData %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanAge = mean(Age),
            seAge = sd((Age)/sqrt(n())))
maturityAgeMean

ggplot(maturityAgeMean, aes(x= Temperature, y= log(meanAge), 
                            group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanAge - seAge, ymax = meanAge + seAge), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Age at Maturity (days)")

# model agemat
ageAccute <- filter(maturityData, Experiment == 'Acute')
ageAcclim <- filter(maturityData, Experiment == "Acclim")

ageMod_accute <- lmer(log(Age) ~ Temperature * Treatment + (1|Clone), data = ageAccute)
ageMod_acclim <- lmer(log(Age) ~ Temperature * Treatment + (1|Clone), data = ageAcclim)
Anova(ageMod_accute, test = "F")
Anova(ageMod_acclim, test = "F")

###### Brood Data Analysis #########
####################################
broodData <- data1 %>% filter(Stage %in% c("B1", "B2", "B3"))
head(broodData)

broodData$Temperature <- factor(broodData$Temperature)

broodDataMean <- broodData %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanBrood = mean(No.Neonates),
            seBrood= sd((No.Neonates)/sqrt(n())))
broodDataMean

# Plot the data

ggplot(broodDataMean, aes(x= Temperature, y= log(meanBrood), 
                          group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  #geom_errorbar(aes(ymin = meanBrood - seBrood, ymax = meanBrood + seBrood), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Average Clutch Size")

### multiple linear regression
clutchMtModlMLR <- lm(No.Neonates ~ Temperature*Experiment + Clone + Treatment, data = broodData)
clutchMtModlMLR 

#Plot the model
autoplot(clutchMtModlMLR)

anova(clutchMtModlMLR)

summary(clutchMtModlMLR)


