library(tidyverse)

library(devtools)
install_github("mdjbru-R-packages/thermPerf")
library(thermPerf)

# FECUNDITY
# sum of 3-clutches

## Source the Trait making script ----

source("./Scripts/Working/MakeFecData_APB.R")

## Fecundity Analyses ----

# look at the data
glimpse(Fec_scale)
p_Fec # the plot

testDat <- filter(Fec_scale, Treatment == 'Control',
                  Experiment == "Acclim",
                  Clone == "D10_A14")

ggplot(testDat, aes(x = Temperature, y = Fec))+
  geom_point()



#testDat <- filter(testDat, Fec>10)

models = getModelLibrary()
names(models)

useModel = getModelLibrary()[c("briere1", "briere2",
                               "candidate01",
                               "candidate02",
                               "candidate03",
                               "candidate04")]

fits = fitModels(useModel, testDat$Temperature, testDat$Fec)
plot(fits)

modDat <- Fec_scale %>% 
  group_by(Treatment, Experiment) %>% 
  nest()

tpc <- function(df){
  useModel = getModelLibrary()[c("briere1", "briere2",
                                 "candidate01",
                                 "candidate02",
                                 "candidate03",
                                 "candidate04",
                                 "candidate05",
                                 "candidate06")]
  fitModels(useModel, df$Temperature, df$Fec)
}

outMods <- modDat %>% 
  mutate(mods = map(data, tpc))

par(mfrow = c(2,2), mar = c(3,2,1,1))
for(i in 1:4){
  plot(outMods$mods[[i]])
}
