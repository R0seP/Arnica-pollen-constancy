#This script aims to investigate whether the specialization index H2 of the 
#pollen-pollinator networks at different sites is influenced by the number of 
#Arnica stems at that site.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#data
species_metrics <- read.csv("species_metrics.csv", h = T)

#model----
#species level d'
#standardized d' against null models
m_d <-  glmmTMB(Arnica_delta_d ~ log(Stems), data = species_metrics)

summary(m_d)

r.squaredGLMM(m_d)

eff_d <- effect("log(Stems)",m_d, xlevels = 50)  
eff.plot(eff_d, plotdata = T,
         ylab = "Blüthgen's d (specialization index) of Arnica",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


#test if model assumptions are met and test model for fit:
qqnorm(resid(m_d)) #looks good
hist(resid(m_d)) #residual distribution looks approaching normal

residuals_d <- simulateResiduals(fittedModel = m_d)
plot(residuals_d)
testOutliers(residuals_d)
#no significance

#effect sizes----
pred_data_d <- data.frame(Stems = c(10, 100, 500))
predictions_d <- predict(m_d, newdata = pred_data_d, type = "response", re.form = NA, se.fit = T)

pred_data_d$Pred.Arnica.d <- predictions_d$fit
pred_data_d$Pred.Arnica.d_SE <- predictions_d$se.fit
pred_data_d
#For 10 stems, the model predicts a d' of 0.38 (+/- 0.06), for 100 stems 
#it predicts a d' of 0.43 (+/-  0.03), for 500 stems a d' of 0.46 (+/- 0.06)

#compare to model with just d (not standardized)
m_comparsion_d <- glmmTMB(Arnica_d ~ log(Stems), data = species_metrics)
summary(m_comparsion_d)
summary(m_d)

eff_comp_d <- effect("log(Stems)",m_comparsion_d, xlevels = 50)  
eff.plot(eff_comp_d, plotdata = T,
         ylab = "Blüthgen's d (specialization index) of Arnica",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)
#as expected, nearly the same since null model d so tiny that standardized
#and not standardized d nearly the same