#libraries
library(tidyverse)
library(nls.multstart)
library(broom)

## source the Trait making script ----
# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are NOT scaled because 0 matters.
# it will produce a plot
source("./Scripts/Working/MakeAllTraits_APB.R")

# doesn't work.
theta_log <- function(x, r, K, theta){
  r*x*(1-(x/K)^theta)
}

## nls.multstart on Subset --------
d_1 <- filter(Fec_scale, 
              Clone == "C14", 
              Treatment == "Control", 
              Experiment == "Acclim") %>% 
  select(Clone, Treatment, Experiment, Temperature, Fec)

nls(Fec ~ r*Temperature*(1-(Temperature/K)^theta),
    data = d_1,
    start = list(K = 25, r = 0.1, theta = 1))

## run nls_multstart with shotgun approach ---
# THIS WORKS WELL ON THIS EXAMPLE, BUT BORKS ON MANY
# TRY C14, PREDATOR, ACCLIM
fit <- nls_multstart(Fec ~ theta_log(r, K, theta, x = Temperature),
                     data = d_1,
                     iter = 250,
                     start_lower = c(r = 0.01, K = 20, theta = 0.5),
                     start_upper = c(r = 0.25, K = 35, theta = 1.2),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(r = 0, K = 0, theta = 1))

tidy(fit)

# predict and plot
newX <- data.frame(Temperature = seq(from = 13, to = 28, length = 100))
out <- augment(fit, newdata = newX) %>% rename(fit = .fitted)
ggplot(out, aes(x = Temperature, y = fit))+
  geom_line()+
  ylab("Fecundity")+
  geom_point(data = d_1, aes(x = Temperature, y = Fec))
