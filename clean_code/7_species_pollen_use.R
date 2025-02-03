#This script examines whether the different pollinator species use Arnica to 
#different amounts.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
library(gridExtra)
library(car)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
comb_all2 <- read.csv("comb_all2.csv", h = T)

comb_all2$Species <- as.factor(comb_all2$Species)

#data visualization----
ggplot(comb_all, aes(fill=Group, y=Nr_Arnica, x=Group)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Sampling group", y = "Nr of Arnica montana pollen in sample")

ggplot(comb_all, aes(fill=Group, y=P_ASTE.Arnica.montana, x=Group)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Sampling group", y = "% Arnica montana pollen in sample")

#calculate the mean P_ASTE.Arnica.montana for each species
species_means_prop <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE)) %>%
  arrange(desc(mean_value))

#reorder the Species factor levels based on the calculated means
comb_all2$Species <- factor(comb_all2$Species, levels = species_means_prop$Species)

#boxplot of percentage of Arnica carried
bp1 <- ggplot(comb_all2, aes(fill = Group, y = P_ASTE.Arnica.montana, x = Species)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12), # increase legend text size
    legend.title = element_text(size = 14), # increase legend title size
    axis.text = element_text(size = 12), # increase axis text size
    axis.title = element_text(size = 14), # increase axis title size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic") # italicize x-axis labels
  ) +
  scale_fill_manual(values = c("area" = "darkgreen", "flower" = "orange")) +
  labs(x = "Species", y = expression("% " * italic("Arnica montana") * " pollen in sample"))



#calculate the mean Nr_Arnica for each species
species_means_nr <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(Nr_Arnica, na.rm = TRUE)) %>%
  arrange(desc(mean_value))

#reorder the Species factor levels based on the calculated means
comb_all2$Species <- factor(comb_all2$Species, levels = species_means_nr$Species)

#boxplot of number of Arnica pollen carried
bp2 <- ggplot(comb_all2, aes(fill = Group, y = Nr_Arnica, x = Species)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12), # increase legend text size
    legend.title = element_text(size = 14), # increase legend title size
    axis.text = element_text(size = 12), # increase axis text size
    axis.title = element_text(size = 14), # increase axis title size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic") # italicize x-axis labels
  ) +
  scale_fill_manual(values = c("area" = "darkgreen", "flower" = "orange")) +
  labs(x = "Species", y = expression("Nr of " * italic("Arnica montana") * " pollen in sample"))

grid.arrange(bp1, bp2, ncol = 1)

#model negative binomial----
#(m9.3 from model selection)

m_species <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                     data = comb_all2, family = nbinom1,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

summary(m_species)

#r-squared value:
r.squaredGLMM(m_species)


#plot predictions for population size
eff_pollen_prop <- effect(c("Stems"),m_species, xlevels = 50)
eff.plot(eff_pollen_prop, plotdata = T,
         ylab = "Nr of Arnica pollen in sample",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, 
         col.data = 3)

#plot predictions for species
eff_species <- effect("Species",m_species, xlevels = 50)

par(mar = c(8.5, 4, 4, 2) + 0.1)  #increase bottom margin
eff.plot(eff_species, plotdata = T,
         ylab = "Nr of Arnica pollen in sample",
         xaxt = "n",
         ylim.data = T, overlay = F, 
         col.data = 3, 
         las = 3)

species_ordered <- as.list(levels(comb_all2$Species))
axis(1, at = 1:22, labels = FALSE)  #suppress default labels
mtext(text = species_ordered, side = 1, at = 1:22, line = 1, las = 2, 
      cex = 0.8, font = 3) #add labels in italcs and rotated

par(mar = c(5, 4, 4, 2) + 0.1) #reset margins



#test if model assumptions are met and test model for fit:
check_overdispersion(m_species) #no overdispersion
check_collinearity(m_species) #no highly correlated predictors

hist(resid(m_species)) #distribution of residuals high peak and more negative

residuals_species <- simulateResiduals(fittedModel = m_species)
plot(residuals_species)
testOutliers(residuals_species)


#analysis of deviance----
#investigate effect of species overall
anova <- Anova(m_species, type = "III")  # Type III ANOVA table
print(anova)


#effect sizes----
#predict the percentage of Arnica pollen for species in different groups for 10, 100, 500 stems
pred_data_species <- data.frame(Stems = c(rep(10,44),rep(100,44),rep(500,44)), 
                                Species = rep(imp_species$Species,6), 
                                Group = c(rep("area",22), rep("flower",22)),
                                nPoll = rep(mean(comb_all2$nPoll), 132))
predictions_species <- predict(m_species, newdata = pred_data_species, type = "response", re.form = NA, se.fit = T)

pred_data_species$Pred.Arnica <- predictions_species$fit
pred_data_species$Pred.Arnica_SE <- predictions_species$se.fit
pred_data_species