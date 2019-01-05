x <- seq(1,45,0.1)
y <- 50*log(x) - 0.1*x^2
plot(y~x, type = "l")

library(tidyverse)

Fec_scale

play <- filter(Fec_scale, Clone == "C14", Treatment == "Predator", 
               Experiment == "Acclim")

ggplot(play, aes(x = Temperature, y = Fec))+
  geom_point()

mod <- lm(Fec ~ Temperature+log(Temperature)+I(Temperature^2), data = play)
summary(mod)
newX <- data.frame(Temperature = seq(13, 28, 1))
newY <- predict(mod, newdata = newX, type = 'response')
pT <- data.frame(newX, fit_Fec = newY)

ggplot(play, aes(x = Temperature, y = Fec))+
  geom_point()+
  geom_line(data = pT, aes(x = Temperature, y = fit_Fec))

Fec_play <- mutate(Fec_scale, 
                   fullFac = paste(Treatment, Experiment, sep = "_"))
play2 <- filter(Fec_play, fullFac == "Control_Acclim")

ggplot(Fec_play, aes(x = Temperature, y = Fec))+
  geom_point()+
  facet_wrap(~fullFac)+
  geom_smooth(method = lm, formula = y ~ x+log(x)+I(x^2))+
  geom_smooth(method = lm, formula = y ~ poly(x,2), colour = 'red', se=FALSE)


ggplot(Fec_scale, aes(x = Temperature, y = Fec, group = Experiment, colour = Experiment))+
  geom_point()+
  facet_grid(Clone~Treatment)+
  geom_smooth(method = lm, formula = y ~ x+log(x)+I(x^2), 
              linetype = 'dashed', se = FALSE)+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE)
