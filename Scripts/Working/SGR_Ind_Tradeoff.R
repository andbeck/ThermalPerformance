library(tidyverse)
library(splines)
library(lme4)
library(car)
library(MESS)
library(lattice)
library(MCMCglmm)

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")

# isolate predation treatment with induction
pred_only <- filter(SGR_Ind, Treatment == "Predator")

# raw data graphs ----------------

# Growth Rate Increases with Temperature
# Acclimated or Acute

# Fig 1 Growth 0 TPC non-linear ---------------------
ggplot(pred_only, aes(x = Temperature, y = Growth0, 
                      colour = factor(Temperature),
                      group = Clone))+
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se=F,
              colour = "grey70")+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  ylab("Somatic Growth Rate (mm/day)")+
  theme_grey(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  guides(color=guide_legend(title="Temp C˚")) 



# Fig 2 Induction Declines or flat with temperature ---------
# Acclim/Acute alterns means

ggplot(pred_only, aes(x = Temperature, y = maxInduction, 
                      colour = factor(Temperature),
                      group = Clone))+
  geom_point(size=2)+
  geom_smooth(method = lm, se=F, col = 'grey70')+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  labs(x = "Temperature (C˚)", y = "Induction Score")+
  theme_grey(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  guides(color=guide_legend(title="Temp C˚")) 


# Trade-off vs. temperature (Raw Data) ----
# mostly negative, with some clones constrained by not inducing
# hint of acclimation-acute change
ggplot(pred_only, aes(x = Growth0, y = maxInduction))+
  geom_smooth(aes(x = Growth0, y = maxInduction), method = lm, se=F,
              colour = 'grey70')+
  geom_jitter(aes(group = Temperature, colour = factor(Temperature)), 
              size = 2, width = 0.1, height = 1)+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(Experiment~Clone)+
  theme_grey(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45))+
  guides(color=guide_legend(title="Temp C˚")) 


# Summarise the Data and Graphs (Raw Data Trade-offs) ----------------

# gather all SGRs (to mat, to clutch 1, to clutch 2)
pred_only2 <- gather(pred_only, key = G_Stage, value = SGR,
                     - ID, -Treatment, -Temperature, -Clone,
                     -Experiment, -maxInduction)

# means of Ind and SGR
mean_SGR_Ind <- pred_only2 %>% 
  group_by(Clone, G_Stage, Temperature, Experiment) %>% 
  summarise(meanInd = mean(maxInduction, na.rm = TRUE),
            meanSGR = mean(SGR, na.rm = TRUE))

# plot all mean Inds ~ Growths (Raw Data Trade-offs) -----
G0 <- filter(mean_SGR_Ind, G_Stage == "Growth0") %>%
  ungroup() %>% 
  mutate(Temperature = factor(Temperature))

# Panel by Experiment, regression among clones
# This shows the overall reduction in SGR with higher Induction, 
# but for acute this is true within ANY temperature
# but for acclimatied, the within temperature is neutral, but the
# trade-offs is negative.
ggplot(G0, aes(x = meanInd, y = meanSGR))+
  # coloured points by temperature
  geom_point(aes(group = Temperature, colour = Temperature))+
  # temperature regressions among clones
  geom_smooth(aes(group = Temperature, colour = Temperature), 
              method = lm, se=FALSE, alpha = 0.5)+
  # overall regression
  geom_smooth(method = lm, se = FALSE, col = 'black')+
  # customise colour scheme
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(Experiment ~ .)

# Panels by Temperature, regression among clones
# This one shows that within any temperature, there is a trade-off
# that varies
# add scales = "free_x" to facet_grid?
ggplot(G0, aes(x = meanInd, y = meanSGR, 
               group = Temperature, colour = Temperature))+
  geom_point(size = 1)+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  geom_smooth(method = lm, se=FALSE, size = 1)+
  facet_grid(Experiment ~ Temperature)+ 
  guides(color=guide_legend(title="Temp C˚"))

# Models of clone means
modTO <- lm(meanSGR ~ meanInd * Temperature * Experiment, data = G0)

Anova(modTO)
summary(modTO)

newX <- expand.grid(
  meanInd = seq(from = 0, to = 100, by = 10),
  Temperature = unique(G0$Temperature),
  Experiment =  unique(G0$Experiment))

newY <- predict(modTO, newdata = newX, interval = 'confidence')

plotThese <- data.frame(newX, newY) %>% 
  rename(predictedSGR = fit, Induction = meanInd)

ggplot(plotThese, aes(x = Induction, y = predictedSGR, 
                      group = Temperature, colour = Temperature))+
  geom_line()+
  geom_point(data = G0, 
             aes(x = meanInd, y = meanSGR))+
  facet_wrap(~Experiment)

# Models: LMER and MGMGglmm ---------------------------------------------

# scale the data
scaleDat <- pred_only %>% as.data.frame() %>% 
  mutate(maxInduction = scale(maxInduction),
         Growth0 = scale(Growth0),
         Growth1 = scale(Growth1)) %>% 
  na.omit()

glimpse(scaleDat)

ggplot(scaleDat, aes(x = maxInduction, y = Growth0, colour = factor(Temperature)))+
  geom_point()+
  facet_grid(Experiment ~ .)
  

# let the effect of growth vary among clone
# let the effect of experiment vary among clone
# let effect of temperature vary among clone

# lmer model ---------------------------------------------------------
modTradeOff <- lmer(Growth0 ~ maxInduction*Temperature*Experiment+
                       (maxInduction|Clone), data = pred_only,
                    control = lmerControl(optimizer = "Nelder_Mead"))

Anova(modTradeOff, test ="F")

# prep plot lmer  Result 

newX <- expand.grid(
  maxInduction = seq(from = 0,
                to = 100,
                length = 10),
  Temperature = unique(pred_only$Temperature),
  Experiment = unique(pred_only$Experiment),
  Clone = unique(pred_only$Clone))

fixed_pred <- predict(modTradeOff, newdata = newX, re.form = NA)
clone_pred <- predict(modTradeOff, newdata = newX, 
                      re.form = ~(1|Clone))

pd <- data.frame(newX, fixed_pred, clone_pred) %>% 
  mutate(Experiment = factor(Experiment, levels = c("Acclim","Acute")))

# lmer model Plot 
lmerFixed <- ggplot(pd, aes(x = maxInduction, y = fixed_pred, colour = factor(Temperature), 
               group = Temperature))+
  geom_line(size = 2)+
  geom_point(data = pred_only, aes(x = maxInduction, y = Growth0), 
             colour = "grey")+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  #scale_y_continuous(breaks = seq(from = -2, to = 2, by = 0.5))+
  facet_grid(Experiment ~ Temperature)+
  labs(y = "SGR", x = "IND")+
  theme_bw(base_size = 15)+
  guides(color=guide_legend(title="Temp C˚"))+
  theme(legend.position = "top")


lmerFixed

lmerClone <- ggplot(pd, aes(x = maxInduction, y = clone_pred, 
                            colour = factor(Temperature), 
                            group = Temperature))+
  geom_line(size = 2)+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(factor(Experiment, levels = c("Acclim","Acute")) ~ Clone)+
  labs(y = "SGR", x = "Ind")+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

lmerClone

# master?
ggplot(pd, aes(x = maxInduction, y = fixed_pred, 
               colour = factor(Temperature), 
               group = Temperature))+
  geom_line(size = 2)+
  scale_colour_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(.~ Experiment)+
  labs(y = "SGR", x = "Ind")+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")


gridExtra::grid.arrange(lmerFixed, lmerClone, ncol = 1)

# MCMCglmm model --------------------------------
# 4 terms in the random effect (intercept + 3)
prior <- list(R = list(V = 1, n = 0.002),
              G = list(G1 = list(V = diag(4), n = 4)))

modBay <- MCMCglmm(maxInduction ~ Growth0*Temperature*Experiment,
                   random = ~us(1+Growth0+Temperature+Experiment):Clone, 
                   family = 'gaussian', prior = prior, pr = TRUE,
                   data = scaleDat)
summary(modBay)$sol

# check fits vs. raw
xyplot(maxInduction + predict(modBay, marginal = NULL) ~
         Growth0|Temperature+Experiment, data = scaleDat)

# create tidy preds
df_plot_bay <- data.frame(scaleDat, 
                          marg_pred = predict(modBay, marginal = NULL),
                          fixed_pred = predict(modBay))

# mcmc Model Plot
mcmcFixed <- ggplot(df_plot_bay, aes(x = Growth0, y = fixed_pred, 
                            colour = Temperature, 
                            group = Temperature))+
  geom_line(size = 2)+
  scale_colour_gradientn(colours = rev(heat.colors(5)))+
  facet_grid(~Experiment)+
  labs(y = "Induction")+
  theme_bw(base_size = 20)


# Compare lmer and MCMC model plot
gridExtra::grid.arrange(lmerFixed, mcmcFixed, ncol = 2)


# marginal plots - but not correct using geom_smooth (see above)
#  Temperature Effect
TEff <- ggplot(df_plot_bay, aes(x = Growth0, y = marg_pred, 
                        colour = Temperature,
                        group = Temperature))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', se = FALSE)+
  scale_colour_gradientn(colours = rev(heat.colors(5)))+
  facet_grid(Experiment ~ .)+
  xlab("Somatic Growth Rate")+ylab("Morphology")+
  theme_bw(base_size = 15)

TEff


CEff <- ggplot(df_plot_bay, aes(x = Growth0, y = marg_pred, 
                                colour = Clone,
                                group = Clone))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm', se = FALSE)+
  #scale_colour_gradientn(colours = rev(heat.colors(5)))+
  facet_grid(Experiment ~ .)+
  xlab("Somatic Growth Rate")+ylab("Morphology")+
  theme_bw(base_size = 15)

CEff

gridExtra::grid.arrange(CEff, TEff, ncol = 2)
