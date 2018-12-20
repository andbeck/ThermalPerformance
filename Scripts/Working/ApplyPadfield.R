# APB Hebe Master 20 Dec 2018

#libraries
library(tidyverse)
library(nls.multstart)
library(broom)

# source the Trait making script ----
# this will make three traits: Growth Rate, Fecundity and max Induction
# all traits are scaled.
# it will produce a plot
source("./Scripts/Working/MakeAllTraits_APB.R")

# Lactin Function ----
lactin <- function(p, Tmax, dT, lam, Temperature){
  exp(p*Temperature) - exp(p*Tmax - ((Tmax - Temperature)/dT)) + lam
}

# nls.multstart on Subset --------
d_1 <- filter(Fec_scale, 
              Clone == "C14", 
              Treatment == "Control", 
              Experiment == "Acclim") %>% 
  select(Clone, Treatment, Experiment, Temperature, Fec)

# run nls_multstart with shotgun approach ---
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
newX <- data.frame(Temperature = seq(from = 13, to = 28, length = 100))
out <- augment(fit, newdata = newX)
ggplot(out, aes(x = Temperature, y = .fitted))+
  geom_line()+
  geom_point(data = d_1, aes(x = Temperature, y = Fec))


# nls.multfit on fuller ----------