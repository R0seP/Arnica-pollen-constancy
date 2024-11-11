#This script aims to investigate whether the specialization index H2 of the 
#pollen-pollinator networks at different sites is influenced by the number of 
#Arnica stems at that site.

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")
library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/EffPlots.R")

#data
network_metrics <- read.csv("network_metrics.csv", h = T)
species_metrics <- read.csv("species_metrics.csv", h = T)

#network level H2----
###model selection----
m1 <- glmmTMB(H2 ~ Stems, data = network_metrics)
m2 <- glmmTMB(H2 ~ log(Stems), data = network_metrics)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m2 with log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m2 with log much higher r2 (though still very low!)

###model----
m_H2 <-  glmmTMB(H2 ~ log(Stems), data = network_metrics)
summary(m_H2)

eff_H2 <- effect("log(Stems)",m_H2, xlevels = 50)  
eff.plot(eff_H2, plotdata = T,
         ylab = "Specialization index H2",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_H2)) #looks nice
hist(resid(m_H2)) #residual distribution looks approximately normal

residuals_H2 <- simulateResiduals(fittedModel = m_H2)
plot(residuals_H2)
testOutliers(residuals_H2)
#no significance


#species level PDI----
###model selection----
m1 <- glmmTMB(Arnica_PDI ~ Stems, data = species_metrics)
m2 <- glmmTMB(Arnica_PDI ~ log(Stems), data = species_metrics)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m1 without log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m1 without log higher r2 (though still very low!)

###model----
m_PDI <-  glmmTMB(Arnica_PDI ~ Stems, data = species_metrics)
summary(m_PDI)

eff_PDI <- effect("Stems",m_PDI, xlevels = 50)  
eff.plot(eff_PDI, plotdata = T,
         ylab = "Paired Differences Index (PDI)",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_PDI)) #some skew
hist(resid(m_PDI)) #residual distribution does not look normal

residuals_PDI <- simulateResiduals(fittedModel = m_PDI)
plot(residuals_PDI)
testOutliers(residuals_PDI)
#no significance


#species level PSI----
###model selection----
m1 <- glmmTMB(Arnica_PSI ~ Stems, data = species_metrics)
m2 <- glmmTMB(Arnica_PSI ~ log(Stems), data = species_metrics)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m2 with log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m2 with log slightly higher r2 (though still very low!)

###model----
m_PSI <-  glmmTMB(Arnica_PSI ~ log(Stems), data = species_metrics)
summary(m_PSI)

eff_PSI <- effect("log(Stems)",m_PSI, xlevels = 50)  
eff.plot(eff_PSI, plotdata = T,
         ylab = "Pollination Service Index (PSI)",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_PSI)) #looks good
hist(resid(m_PSI)) #residual distribution does not look normal
#(much bigger residuals at higher x values)

residuals_PSI <- simulateResiduals(fittedModel = m_PSI)
plot(residuals_PSI)
testOutliers(residuals_PSI)
#no significance


#species level d'----
#standardized d' against null models
###model selection----
m1 <- glmmTMB(Arnica_delta_d ~ Stems, data = species_metrics)
m2 <- glmmTMB(Arnica_delta_d ~ log(Stems), data = species_metrics)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m2 with log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m2 with log slightly higher r2 (though still very low!)

###model----
m_d <-  glmmTMB(Arnica_delta_d ~ log(Stems), data = species_metrics)
summary(m_d)

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