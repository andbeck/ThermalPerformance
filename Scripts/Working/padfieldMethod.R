library(tidyverse)
library(splines)
library(mgcv)
library(broom)

# data
hd <- read_csv("./Data/PopulationGrowth.csv")
glimpse(hd)


# Get Ro data in order
# including standardisation
Ro <- hd %>%
  filter(Stage == "B1"|Stage == "B2") %>% 
  mutate(lxmx = Survival * Mean.Rep) %>% 
  ungroup() %>% 
  group_by(Clone, Temperature, Treatment, Experiment) %>% 
  summarise(Ro = sum(lxmx, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # create combined factor, unordered
  mutate(TrExp = factor(paste(Treatment,Experiment, sep = "_"), ordered = FALSE)) %>% 
  # ensure clone is unordered factor (from paper)
  mutate(Clone = factor(Clone, ordered = FALSE)) %>% 
  # ensure Temperature is numeric
  mutate(Temperature = as.numeric(Temperature)) %>% 
  mutate(Ro = scale(Ro))

# first view of it
ggplot(Ro, aes(x = Temperature, y = Ro, group = Clone, colour = Clone))+
  geom_point()+
  geom_smooth(span = 2, se=FALSE)+
  facet_grid(Treatment ~ Experiment)

# padfield ----------
devtools::install_github("padpadpadpad/nls.multstart")
library(nls.multstart)

# function
# apb lactin from Luhrling
lactin <- function(p, Tmax, dT, lam, Temperature){
  exp(p*Temperature) - exp(p*Tmax - ((Tmax - Temperature)/dT)) + lam
}

d_1 <- filter(Ro, Treatment == "Control" & Experiment == "Acute")

ggplot(d_1, aes(x = Temperature, y = Ro))+
         geom_point()+
         geom_smooth(method = lm, formula = y ~ poly(x,2), se = FALSE)

# run nls_multstart with shotgun approach
fit <- nls_multstart(Ro ~ lactin(p, Tmax, lam, dT, Temperature = Temperature),
                     data = d_1,
                     iter = 250,
                     start_lower = c(p = -1, Tmax = 10, dT = 1, lam = -2),
                     start_upper = c(p = 1, Tmax = 30, dT = 10, lam = 1),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(p = -10, Tmax = 0, dT = 0, lam = 0))

tidy(fit)
newX <- data.frame(Temperature = seq(13, 28, 1))
out <- augment(fit, newdata = newX)
ggplot(out, aes(x = Temperature, y = .fitted))+
  geom_line()+
  geom_point(data = d_1, aes(x = Temperature, y = Ro))

