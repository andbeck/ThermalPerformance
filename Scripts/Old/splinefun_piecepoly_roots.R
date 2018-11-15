# EXAMPLE OF FITTING SPLINE TO PREDICTION TO GET VERTICES ------------------------
# with Mixed Model

# test one of the clones predictions
grabOne <- predDat %>% filter(Clone == "LD33" & Treatment == "Control" &
                                Experiment == "Acute")

# convert predictions to spline mod
splineMod <- with(grabOne, smooth.spline(pred_fix ~ Temperature))

# calculate derivative to get turning points
p_deriv <- data.frame(predict(splineMod, deriv = 1))

# grab the temperature (x) at which turn happens
# this finds the derivative value closest to 0
Topt <- p_deriv[which.min(abs(p_deriv$y - 0)),]$x

# grab the trait - value (y) at which turn happens
Popt <- predict(splineMod, x = Topt)$y

# plot the subset
ggplot(grabOne, aes(x = Temperature, y = pred_fix))+
  geom_line()+
  geom_vline(aes(xintercept = Topt), linetype = "dashed")+
  geom_hline(aes(yintercept = Popt), linetype = "dashed")

# create function to get Topt and Popt -----------------------------------------------------------
# note addition to deal with clones that don't have turning point
# must find + and - derivative values (gradients) to have turning point
# TODO ADD AUC

TPOpt <- function(predData){
  splineMod <- with(predData, smooth.spline(pred_fix ~ Temperature))
  p_deriv <- data.frame(predict(splineMod, deriv = 1))
  
  # test for turning point
  if(all(p_deriv>0)) {
    Topt <- last(p_deriv$x)
    Popt <- last(predict(splineMod)$y)
  } else {
    Topt <- p_deriv[which.min(abs(p_deriv$y - 0)),]$x
    Popt <- predict(splineMod, x = Topt)$y  
  }
  return(data.frame(Topt, Popt))
}

# USE function against all predictions to get Opts -----------------------------

# creat group structure of 40 sets from PREDICTIONS of lmer model
byPEC <- predDat %>% 
  group_by(Treatment, Experiment, Clone) %>% 
  nest()

# get all the models using mutate(map) tricks
outMods <- byPEC %>% 
  mutate(Opts = map(data, TPOpt))

# unnest the tibble and gather the stuff into approporiate columns
outDats <- outMods %>% unnest(Opts) %>% 
  select(-data) %>% 
  gather(OptCode, OptVal, -Treatment, -Experiment, -Clone) %>% 
  data.frame()

# summarise these derived values among clones
sumPlot <- outDats %>% 
  group_by(Experiment, Treatment, OptCode) %>% 
  summarise(meanOptVal = mean(OptVal),
            seOptVal = sd(OptVal)/sqrt(10))

# plot the derived values
ggplot(sumPlot, aes(x = Experiment, y = meanOptVal, 
                    colour = Treatment, group = Treatment))+
  geom_line(position = position_dodge(0.1))+
  geom_errorbar(aes(ymin = meanOptVal - seOptVal, ymax = meanOptVal + seOptVal),
                width = 0.1,
                position = position_dodge(0.1))+
  facet_wrap(~OptCode, scales = 'free_y')

# fit models to derived values
# effect size and anova calculation valuable
outDats %>% 
  filter(OptCode == "Popt") %>% 
  lm(OptVal ~ Treatment * Experiment, data = .) %>% 
  etasq(., anova = TRUE)

outDats %>% 
  filter(OptCode == "Topt") %>% 
  lm(OptVal ~ Treatment * Experiment, data = .) %>% 
  etasq(., anova = TRUE)


# b-spline with PiecePoly - not possible with mixed model -----------
mod_bs <- lm(Growth0 ~ bs(Temperature), data = ConAcute)

getRoots <- function(splineMod){
  require(SplinesUtils, quietly = TRUE)
  mm <- RegBsplineAsPiecePoly(splineMod, SplineTerm = "bs(Temperature)")
  Topt <- solve(mm, deriv = 1)
  Popt <- predict(mod_bs, newdata = data.frame(Temperature = Topt))
  return(list(Topt = Topt, Popt = Popt))
}

rr <- getRoots(mod_bs)

# plot and predictions from models
plot(Growth0 ~ Temperature, data = ConAcute)
preds <- predict(mod_bs, newdata = newX)
lines(newX$Temperature, preds, col = 4)
abline(v = rr$Topt, h = rr$Popt, lty = 3)

