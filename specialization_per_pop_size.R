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

#model selection----
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

#model
m_specialization <-  glmmTMB(H2 ~ log(Stems), data = network_metrics)
summary(m_specialization)

eff_specialization <- effect("log(Stems)",m_specialization, xlevels = 50)  
eff.plot(eff_specialization, plotdata = T,
         ylab = "Specialization index H2",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)
