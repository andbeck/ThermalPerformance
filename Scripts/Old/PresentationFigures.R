library(tidyverse)

# data
SGR_Ind <- read_csv("./Data/growthInd_HC0.csv")
Ro <- read_csv("./Data/Ro_Data.csv")

# isolate predation treatment with induction
pred_only <- filter(SGR_Ind, Treatment == "Predator")
pred_only_Ro <- filter(Ro, Treatment == "Predator")

meanInd <- pred_only %>% 
  group_by(Clone, Temperature, Experiment) %>% 
  summarise(meanInd = mean(maxInduction, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(Clone, Temperature, Experiment)

meanRo <- pred_only %>% select(Clone, Temperature, Experiment, Growth0) %>%
  group_by(Clone, Temperature, Experiment) %>% 
  summarise(meanRo = mean(Growth0, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(Clone, Temperature, Experiment)

head(meanRo)
head(meanInd)
Ro_Ind_Trade <- left_join(meanRo, meanInd)


# Clone x Treatment x Experiment Plot ----

# SGR ----
ggplot(SGR_Ind, aes(x = Temperature, y = Growth1, colour = Treatment,
                    group = Treatment))+
  geom_point(size = 2)+
  geom_smooth(span = 1.5, se=F)+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  ylab("Somatic Growth Rate")+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45))

# SGR - C14 single clone
C14 <- filter(SGR_Ind, Clone == "C14")

ggplot(C14, aes(x = Temperature, y = Growth1, 
                colour = Treatment, group = Treatment))+
  geom_point(size = 2)+
  geom_smooth(span = 1.5, se=F)+
  facet_grid(Experiment~.)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  ylab("Somatic Growth Rate")+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45))

# PGR ----
ggplot(Ro, aes(x = Temperature, y = Ro.B2, colour = Treatment,
               group = Treatment))+
  geom_point(size = 2) + ylim(0,50) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x,2))+
  facet_grid(Experiment ~ Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  ylab("Population Growth Rate")+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45))

# Induction ----
# reorder Clone by mean of maxInduction
ggplot(pred_only, aes(x = reorder(Clone, maxInduction, median), 
                      y = maxInduction, fill = Clone))+
  geom_boxplot()+
  facet_grid(Experiment ~ .)+
  ylab("Induction Score")+
  xlab("Clone")+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
        legend.position = "none")

# Induction Declines or flat with temperature
# Acclim/Acute alterns means
ggplot(pred_only, aes(x = Temperature, y = maxInduction, group = Clone,
                      colour = Clone))+
  geom_jitter(size=2, height = 0, alpha = 0.6)+
  geom_smooth(method = lm, se=F)+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
        legend.position = "none")


ggplot(SGR_Ind, aes(x = Temperature, y = maxInduction, group = Treatment,
                      colour = Treatment))+
  geom_jitter(size=2, height = 0, alpha = 0.6)+
  geom_smooth(method = lm, se=F)+
  facet_grid(Experiment~Clone)+
  scale_x_continuous(breaks = c(13,16,20,24,28))+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))


## Ro vs Induction (real cost?) ----
ggplot(Ro_Ind_Trade, aes(x = meanInd, y = meanRo, 
                         colour = Temperature, 
                         group = Experiment))+
  geom_point(size = 3)+
  scale_colour_gradientn(colours = rev(heat.colors(5)))+
  geom_smooth(span = 2, se = FALSE)+
  xlab("Morphology Defence") + ylab("Population Growth Rate")+
  facet_grid(Experiment~.)+
  theme_bw(base_size = 15)

# + ylim(0,50) +
#   geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x,2))+
#   facet_grid(Experiment ~ Clone)+
#   scale_x_continuous(breaks = c(13,16,20,24,28))+
#   ylab("Population Growth Rate")+
#   theme_bw(base_size = 15)+
#   theme(axis.text.x = element_text(angle = 45))
