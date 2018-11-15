library(tidyverse)
library(splines)
library(nlme) # needed for nlsList
library(lme4)
library(car)
library(pastecs)
#devtools::install_github('ZheyuanLi/SplinesUtils')
library(SplinesUtils)

# TODO :::: diagnose clones that mean something and those that fuck up/
# TODO :::: figure out nlme mixed?
# TODO :::: stick with fixed effects model

# data
Ro <- read_csv("./Data/Ro_Data.csv")
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")

## SGR ------------------------------------------

# create Factorial column that is 
# Treatment:Experiment:Clone grouping
# Create a scaled dataset and a response one with Factorial

scaleDat <- SGR_Ind %>% as.data.frame() %>% 
  mutate(maxInduction = scale(maxInduction),
         Growth0 = scale(Growth0),
         Growth1 = scale(Growth1),
         factorial = paste(Treatment, Experiment, Clone, sep = ":")) %>%
  na.omit() %>% 
  data.frame()

SGR_Ind2 <- SGR_Ind %>% 
  mutate(factorial = paste(Treatment, Experiment, Clone, sep = ":")) %>%
  na.omit()

# first plot of raw or scaled data
ggplot(SGR_Ind, aes(x = Temperature, y = Growth0, 
                                   colour = Treatment))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)+
  facet_grid(Experiment ~ Clone)+
  theme_bw()


# can try and fit B-splines as mixed effect and convert using 
# SplineUtils or fit with nlme/nls vertex model

# test nls on a subset ------------------------------------------------------------------------
# Control, Acute and Clone LD33
ConAcute <- subset(SGR_Ind2, Treatment == "Control" & 
                     Experiment == "Acute" & 
                     Clone == "LD33")
ggplot(ConAcute, aes(x= Temperature, y = Growth0))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2))

# vertex model of polynomial
m1 <- nls(Growth0 ~ a*(Temperature - h)^2 + K,
          start = list(a = -1, h = 25, K = 1),
          data = ConAcute)

# nlme random effects method failing
# can't seem to get model to fit any covariates to fixed.

mod <- nlme(Growth0 ~ a*(Temperature - h)^2 + K,
            fixed = a+h+K~1,
            random = h+K ~ 1|Clone,
            start = list(fixed= c(a = -1, h = 25, K = 1)),
            data = gd)

summary(mod)

# This is the vertex formula for the quadratic ------------------------------------------
# h and K are the x and y coordinates of the turning point
# and the lactin1 and lactin2 models.
# in lactin1, p is supplsed to be max performance? (not)
# in lactin2, p and LL determine peak and slope of ascendency

vertex <- function(a, h, K, Temperature){a*(Temperature - h)^2 + K}
lactin1 <- function(p, Tmax, dT, Temperature){
  exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT)}
lactin2 <- function(p, Tmax, dT, LL, Temperature){
  exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT) + LL}

# for example, this vertex peaks at x = 25 and a value of y = 1
# numerical solutions required for the lactins

TT <- seq(13,30,1)
VV_p <- vertex(a = -0.01, h = 25, K = 1, TT)
VV_lac1 <- lactin1(p = 0.1, Tmax = 28, dT = 6, Temperature = TT)
VV_lac2 <- lactin2(p = 0.1, Tmax = 28, dT = 6, LL = -1, Temperature = TT)

df <- data.frame(TT,VV_p, VV_lac1, VV_lac2)
ggplot(df, aes(x = TT, y = VV_p))+
  geom_line()+
  geom_line(aes(y = VV_lac1), col = 'red')+
  geom_line(aes(y = VV_lac2), col = 'green')+
  ylim(0,4)

# test nls on a subset ------------------------------------------------------------------------
# Control, Acute and Clone LD33
ConAcute <- subset(SGR_Ind2, Treatment == "Control" & 
                     Experiment == "Acute" & 
                     Clone == "LD33")
ggplot(ConAcute, aes(x= Temperature, y = Growth0))+
  geom_point()+
  geom_smooth(method = lm, formula = y ~ poly(x,2))

# vertex model of polynomial
m1 <- nls(Growth0 ~ a*(Temperature - h)^2 + K,
    start = list(a = -1, h = 25, K = 1),
    data = ConAcute)

# lactin1 model
m2 <- nls(Growth0 ~ exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT),
          start = list(p = 0.1, Tmax = 28, dT = 6),
          data = ConAcute)

# lactin2 model
m3 <- nls(Growth0 ~ exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT) + LL,
          start = list(p = 0.1, Tmax = 28, dT = 2, LL = -2),
          data = ConAcute)

m4 <- lm(Growth0 ~ bs(Temperature), data = ConAcute)

# look at coefficients
coef(m1); coef(m2); coef(m3); coef(m4)

# plot and predictions from models
plot(Growth0 ~ Temperature, data = ConAcute)
newX <- data.frame(Temperature = seq(13,28, length = 50))
preds1 <- predict(m1, newdata = newX)
preds2 <- predict(m2, newdata = newX)
preds3 <- predict(m3, newdata = newX)
preds4 <- predict(m4, newdata = newX)
lines(newX$Temperature, preds1, col = 1)
lines(newX$Temperature, preds2, col = 2)
lines(newX$Temperature, preds3, col = 3)
lines(newX$Temperature, preds4, col = 4)

# how to find derivative using splines -------------------------------------
# see B-SplinesMethod Script for full use 
# via Purrr.
devtools::install_github('ZheyuanLi/SplinesUtils')
library(SplinesUtils)

mm <- RegBsplineAsPiecePoly(m4, SplineTerm = "bs(Temperature)")
solve(mm, deriv = 1)

# use nlsList for multiple fits ----------------------------------------

# calculate 40 estimates (2 x 2 x 10 clones) 
# 2 treatments, 2 experiment exposures and 10 clones

outVertex <- nlsList(Growth0 ~ a*(Temperature - h)^2 + K|factorial,
        start = list(a = -1, h = 25, K = 2),
        data = SGR_Ind2)

outOffset <- lmList(Growth0 ~ poly(Temperature, 3)|factorial, scaleDat)

ggplot(SGR_Ind2, aes(x = Temperature, y = Growth0))+
  geom_point()+
  facet_wrap(~factorial, ncol = 10)

# diagnose fuck ups
ConAccCy <- subset(SGR_Ind2, Treatment == "Control" & 
                     Experiment == "Acclim" & 
                     Clone == "Cyril")

# vertex model ConAccCy
mConAccCy <- nls(Growth0 ~ a*(Temperature - h)^2 + K,
                 start = list(a = -1, h = 25, K = 1),
                 data = ConAccCy)

coef(mConAccCy) # temperature of 0.8 and peak of 0.003 ???

# failure clear
plot(Growth0 ~ Temperature, data = ConAccCy)
newX <- data.frame(Temperature = seq(13,28, length = 50))
preds1 <- predict(m1, newdata = newX)
lines(preds1)


# PURRR approach to multiple models (lmList like) ------------------------------
# created group structure of 40 sets
byPEC <- SGR_Ind %>% 
  group_by(Treatment, Experiment, Clone) %>% 
  nest()

# generic model
polyMod <- function(df){
  nls(Growth0 ~ a*(Temperature - h)^2 + K,
      start = list(a = -1, h = 25, K = 1), 
      data = df)}

# get all the models
outMods <- byPEC %>% 
  mutate(Models = map(data, polyMod),
         sumTables = Models %>% map(tidy),
         coefs = sumTables %>% map("estimate"),
         Topt = coefs %>% map_dbl(2),
         Popt = coefs %>% map_dbl(3))

outMods

ggplot(outMods, aes(x = Experiment, y = Topt, colour = Treatment))+
  geom_point()+
  facet_wrap(~Treatment)


# lactin1 approach fails even more ---------------------------------------
outlactin1 <- nlsList(Growth0 ~ exp(p*Temperature) - exp(p*Tmax - (Tmax - Temperature)/dT)|factorial,
                      start = list(p = 0.1, Tmax = 28, dT = 6),
                      data = SGR_Ind2)

# manipulate the lmList coefficients into a data frame  ----------------------------------
# separate the factorial back into the appropriate things.
df <- data.frame(coef(outDat)) %>% 
  rownames_to_column(var = "ID") %>% 
  separate("ID", into = c("Treatment", "Experiment", "Clone") )

# create a summary version for plotting
dfSum <- df %>% 
  group_by(Treatment, Experiment) %>% 
  summarise(
    aM = mean(a, na.rm = TRUE),
    aS = sd(a, na.rm = TRUE)/sqrt(sum(!is.na(a))),
    hM = mean(h, na.rm = TRUE),
    hS = sd(h, na.rm = TRUE)/sqrt(sum(!is.na(h))),
    kM = mean(K, na.rm = TRUE),
    kS = sd(K, na.rm = TRUE)/sqrt(sum(!is.na(K)))
  )

AA <- ggplot(dfSum, aes(x = Experiment, y = aM, colour = Treatment,
                  group = Treatment,
                  ymin = aM - aS, ymax = aM + aS))+
  geom_point(position = position_dodge(0.1))+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(width = 0.1, position = position_dodge(0.1))

HH <- ggplot(dfSum, aes(x = Experiment, y = hM, colour = Treatment,
                        group = Treatment,
                        ymin = hM - hS, ymax = hM + hS))+
  geom_point(position = position_dodge(0.1))+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(width = 0.1, position = position_dodge(0.1))

KK <- ggplot(dfSum, aes(x = Experiment, y = kM, colour = Treatment,
                        group = Treatment,
                        ymin = kM - kS, ymax = kM + kS))+
  geom_point(position = position_dodge(0.1))+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(width = 0.1, position = position_dodge(0.1))

# visualise 
# NOT MUCH GOING ON?
gridExtra::grid.arrange(HH,KK)

# manova of x and y
# NOWT
mod <- lm(cbind(h,K) ~ Treatment * Experiment, data = na.omit(df))
Anova(mod)

# anovas of x and y
# NOWT

mod <- lm(h ~ Treatment * Experiment, data = na.omit(df))
Anova(mod)




