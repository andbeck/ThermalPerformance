# Hebe All Data
rm(list=ls())
library(tidyverse)

# all data
#allDat <- read_csv("AllResults_17_07_18.csv")
allDat <- read_csv(file.choose())

# check your data and levels of the factors
distinct(allDat, Clone)
distinct(allDat, Treatment)
distinct(allDat, Temperature)
distinct(allDat, Experiment)
distinct(allDat, Stage) # problems
distinct(allDat, Mature) # problems

# fix the factor levels
allDat2 <- 
  allDat %>%
  mutate(Experiment = ifelse(Experiment == "Aclim", "Acclim", Experiment)) %>% 
  mutate(Stage = ifelse(Stage == "Mature", "mt", Stage)) %>% 
  mutate(Clone = ifelse(Clone == "LD34"| Clone == "LD35", "LD33", Clone))

# check
distinct(allDat2, Stage)
distinct(allDat2, Clone)

# get the mature data
matDat <- filter(allDat2, Stage == "mt")
matDat
# plot it to isolate other issues
ggplot(matDat, aes(x = factor(Temperature), y = Body, 
                   colour = Treatment, group = Treatment))+
  geom_jitter(height = 0, width = 0.25)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()

# consider removing D10_74
# D10_45 has some very low numbers at 13C


# check sample sizes by temp and treatment
# something seems fishy here.
sampSize <- matDat %>% 
  group_by(Experiment, Treatment, Temperature, Clone) %>% 
  #summarise(nums = sum(!is.na(Body)))
  summarise(nums = n())


ggplot(sampSize, aes(x = factor(Temperature), y = nums, 
                     colour = Treatment, group = Treatment))+
  geom_point()+
  geom_line()+
  facet_grid(Experiment ~ Clone)+
  scale_y_continuous(breaks = 1:10)+
  theme_bw()


# Age Plot - isloate other issues
ggplot(matDat, aes(x = factor(Temperature), y = Age, 
                   colour = Treatment, group = Treatment))+
  geom_jitter(height = 0, width = 0.25)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()


#### Calculate the means ####
#Size
matDatMean <- matDat %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>% 
  summarise(meanSize = mean(Body, na.rm = TRUE),
            seSize = sd(Body, na.rm = TRUE)/
              sqrt(sum(!is.na(Body))))

# seperated clones
ggplot(matDatMean, aes(x= Temperature, y= meanSize, group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanSize - seSize, ymax = meanSize + seSize), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Size at Maturity (mm)")

# all clones together 
ggplot(matDatMean, aes(x= Temperature, y= meanSize, group = Clone, colour = Clone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanSize - seSize, ymax = meanSize + seSize), width = .1) +
  facet_grid(Experiment ~ Treatment) +
  ylab("Mean size at maturity (mm)")

### Age
matAgeMean <- matDat %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanAge = mean(Age, na.rm=TRUE),
            seAge = sd(Age, na.rm = TRUE)/sqrt(n()))

# plot the data - seperated clones
ggplot(matAgeMean, aes(x= Temperature, y= meanAge, group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanAge - seAge, ymax = meanAge + seAge), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Age at Maturity (mm)")

# All clones together
ggplot(matAgeMean, aes(x= Temperature, y= meanAge, group = Clone, colour = Clone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanAge - seAge, ymax = meanAge + seAge), width = .1) +
  facet_grid(Experiment ~ Treatment) +
  ylab("Age at Maturity (mm)")

###### Clutch Size ####
#Average per clutch:
broodDat <- allDat %>% filter(Stage %in% c("B1", "B2", "B3"))
broodDat

# look for outliers
ggplot(broodDat, aes(x = factor(Temperature), y = No.Neonates, 
                   colour = Treatment, group = Treatment))+
  geom_jitter(height = 0, width = 0.25)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()


### Average clutch size #####
broodDatMean <- broodDat %>%
  group_by(ID, Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanBrood = mean(No.Neonates, na.rm = TRUE),
            seBrood= sd(No.Neonates, na.rm = TRUE)/sqrt(n()))

broodDatMean2 <- broodDatMean %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanBrood2 = mean(meanBrood, na.rm = TRUE),
            seBrood2= sd(meanBrood, na.rm = TRUE)/sqrt(n()))

#Plot the data - clones seperated
ggplot(broodDatMean2, aes(x= Temperature, y= meanBrood2, group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanBrood2 - seBrood2, ymax = meanBrood2 + seBrood2), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Average clutch size")

# all clones together
ggplot(broodDatMean2, aes(x= Temperature, y= meanBrood2, group = Clone, colour = Clone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanBrood2 - seBrood2, ymax = meanBrood2 + seBrood2), width = .1) +
  facet_grid(Experiment ~ Treatment) +
  ylab("Average clutch size")

### sum of clutch size ###
broodDatSum <- broodDat %>%
  group_by(ID, Clone, Experiment, Treatment, Temperature) %>%
  summarise(sumBrood = sum(No.Neonates, na.rm = TRUE),
            seBrood= sd(sumBrood, na.rm = TRUE)/sqrt(n()))


broodDatMeanSum <- broodDatSum %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>%
  summarise(meanSumBrood = mean(sumBrood, na.rm = TRUE),
            seSumBrood= sd(sumBrood, na.rm = TRUE)/sqrt(n()))


ggplot(broodDatMeanSum, aes(x= Temperature, y= meanSumBrood, group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanSumBrood - seSumBrood, ymax = meanSumBrood + seSumBrood), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Total No. Neonates")

# all clones together
ggplot(broodDatMeanSum, aes(x= Temperature, y= meanSumBrood, group = Clone, colour = Clone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanSumBrood - seSumBrood, ymax = meanSumBrood + seSumBrood), width = .1) +
  facet_grid(Experiment ~ Treatment) +
  ylab("Total No. Neonates")

# Different number of broods between acute and acclim so 
# can't compare (acclim had 3 broods but acute only had 2)
# some may have died and so may only have one clutch


###### Growth Rate ######
GDat <- filter(allDat, Age == "1" | Stage == "mt")
GDat

initialDat <- GDat %>% mutate(ISize = (Age == "1" ))
initialDat

growth <- initialDat %>% group_by(ID) %>% arrange(Age) %>% 
  mutate(Growth = log(last(Body)/ first(Body))/ last(Age))

GrowthRate <- filter(growth, ISize == "FALSE")


# Mean growth rate
GrowthDatMean <- GrowthRate %>%
  group_by(Clone, Treatment, Temperature, Experiment) %>% 
  summarise(meanGrowth = mean(Growth, na.rm = TRUE),
            seGrowth = sd(Growth, na.rm = TRUE)/
              sqrt(sum(!is.na(Growth))))

# clones seperated
ggplot(GrowthDatMean, aes(x= Temperature, y= meanGrowth, group = Treatment, colour = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanGrowth - seGrowth, ymax = meanGrowth + seGrowth), width = .1) +
  facet_grid(Experiment ~ Clone) +
  ylab("Average growth rate (mm per day)")

### All clones together
ggplot(GrowthDatMean, aes(x= Temperature, y= meanGrowth, group = Clone, colour = Clone)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = meanGrowth - seGrowth, ymax = meanGrowth + seGrowth), width = .1) +
  facet_grid(Experiment ~ Treatment) +
  ylab("Average growth rate (mm per day)")
