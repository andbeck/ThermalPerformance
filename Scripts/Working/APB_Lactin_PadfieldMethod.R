# APB Hebe Lactin 20 Dec 2018

#libraries
library(tidyverse)
library(nls.multstart)
library(broom)

## source the Trait making script ----
# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are NOT scaled because 0 matters.
# it will produce a plot
source("./Scripts/Working/MakeAllTraits_APB.R")

## Lactin Function ----
# used by Luhring and deLong

lactin <- function(p, Tmax, dT, lam, Temperature){
  exp(p*Temperature) - exp(p*Tmax - ((Tmax - Temperature)/dT)) + lam
}

## nls.multstart on Subset --------
d_1 <- filter(Fec_scale, 
              Clone == "C14", 
              Treatment == "Control", 
              Experiment == "Acclim") %>% 
  select(Clone, Treatment, Experiment, Temperature, Fec)

## run nls_multstart with shotgun approach ---
# THIS WORKS WELL ON THIS EXAMPLE, BUT BORKS ON MANY
# TRY C14, PREDATOR, ACCLIM
fit <- nls_multstart(Fec ~ lactin(p, Tmax, lam, dT, Temperature = Temperature),
                     data = d_1,
                     iter = 250,
                     start_lower = c(p = -1, Tmax = 10, dT = 1, lam = -2),
                     start_upper = c(p = 1, Tmax = 30, dT = 10, lam = 1),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(p = -10, Tmax = 0, dT = 0, lam = 0))

tidy(fit)

# predict and plot
newX <- data.frame(Temperature = seq(from = 13, to = 28, length = 100))
out <- augment(fit, newdata = newX) %>% rename(fit = .fitted)
ggplot(out, aes(x = Temperature, y = fit))+
  geom_line()+
  ylab("Fecundity")+
  geom_point(data = d_1, aes(x = Temperature, y = Fec))

## get parameters Function
getParams <- function(model, Tmin = -10, Tmax = 35, nit = 1000){
  newX <- data.frame(Temperature = seq(from = Tmin, to = Tmax, length = nit))
  out <- augment(model, newdata = newX) %>% 
    rename(fit = .fitted)
  CTmax_min <- out %>% filter(abs(fit - 1) == min(abs(fit - 1)))
  Pmax <- max(out$fit)
  Topt <- filter(out, fit == as.numeric(Pmax))$Temperature
  return(list(CTmax_min = CTmax_min,
              PerfMax = Pmax,
              Topt = Topt))
  
}

getParams(fit)

# nls.multfit on fuller ----------
