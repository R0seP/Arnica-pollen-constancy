##This script explores the influence of Arnica population size on the pollinator
#community composition


library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
comb1 <- read.csv("comb1.csv", h = T)

###model selection----
m1 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems), 
                     data = comb1, family = binomial)

m2 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ Stems, 
           data = comb1, family = binomial)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m1 with log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m1 with log also higher r2, use that model!


###model----
#model the impact the size of an Arnica patch has on the pollinator community by
#using a binomial model with Arnica-associated species as "success", and area-associated
#species as "failure", without the species "Meligethes":

m_community <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems), 
                       data = comb1, family = binomial)
summary(m_community)

#r-squared
r.squaredGLMM(m_community)

eff_community <- effect("log(Stems)",m_community, xlevels = 50)  
eff.plot(eff_community, plotdata = T,
         ylab = "Proportion of Arnica-associated pollinators",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)



#test if model assumptions are met and test model for fit:
check_overdispersion(m_community)

hist (resid(m_community))

residuals_community <- simulateResiduals(fittedModel = m_community)
plot(residuals_community)
testOutliers(residuals_community)
#residual vs predicted significant deviations


#alternative model----
#try model with quadratic component to improve fit:
#selection: 
m1 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems), 
              data = comb1, family = binomial)
m2 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems)
              + I((Stems)^2), data = comb1, family = binomial)
mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m1 with quadratic term ranked slightly higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#much higher r2 of model with quadratic term

summary(m2)

#plot model with quadratic component
eff_community2_quad <- effect("log(Stems)",m2, xlevels = 50)  
eff.plot(eff_community2_quad, plotdata = T,
         ylab = "Proportion of Arnica-associated pollinators",
         xlab = "Population size Arnica (Nr Stems)",
         main = "community composition without Meligethes",
         ylim.data = T, overlay = F, col.data = 3)
#looks overfitted!


###effect sizes----
#r-squared value:
r.squaredGLMM(m_community2)

pred_data_comm2 <- data.frame(Stems = c(10, 100, 500))
predictions_comm2 <- predict(m_community2, newdata = pred_data_comm2, type = "response", se.fit = TRUE)

pred_data_comm2$Pred.Arnica <- predictions_comm2$fit
pred_data_comm2$Pred.Arnica_SE <- predictions_comm2$se.fit
pred_data_comm2
#for 10 stems, 36.0 (+/- 4.5)% of Arnica pollinators are expected, for 100 stems
#47.6 (+/- 2.5)% of Arnica pollinators, and for 500 stems 55.8 (+/- 4.5)% of 
#Arnica associated pollinators are expected.