##This script explores the influence of Arnica population size on the pollinator
#community composition


setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")
library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/EffPlots.R")

#get data
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)

###model selection----
m1 <- glmmTMB(cbind(Arnica_associated, Area_associated) ~ log(Stems), 
                     data = comb1, family = binomial)

m2 <- glmmTMB(cbind(Arnica_associated, Area_associated) ~ Stems, 
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
#model the impacct the size of an Arnica patch has on the pollinator community by
#using a binomial model with Arnica-associated species as "success", and area-associated
#species as "failure":

m_community <- glmmTMB(cbind(Arnica_associated, Area_associated) ~ log(Stems), 
                       data = comb1, family = binomial)
summary(m_community)
check_overdispersion(m_community)
#no overdispersion, keep current model

eff_community <- effect("log(Stems)",m_community, xlevels = 50)  
eff.plot(eff_community, plotdata = T,
         ylab = "Proportion of Arnica-associated pollinators",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
#qqnorm(resid(m_community)) #qqplot looks good
qqnorm(resid(m_community))

residuals_community <- simulateResiduals(fittedModel = m_community)
plot(residuals_community)
testOutliers(residuals_community)
#no outliers or significant deviations detected
#model seems like a good fit overall!


###effect sizes----
#r-squared value:
r.squaredGLMM(m_community)


#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_comm <- data.frame(Stems = c(10, 100, 500))
predictions_comm <- predict(m_community, newdata = pred_data_comm, type = "response", se.fit = TRUE)

pred_data_comm$Pred.Arnica <- predictions_comm$fit
pred_data_comm$Pred.Arnica_SE <- predictions_comm$se.fit
pred_data_comm
#for 10 stems, 44.2 (+/- 4.5)% of Arnica pollinators are expected, for 100 stems
#51.3 (+/- 2.4)% of Arnica pollinators, and for 500 stems 56.3 (+/- 4.3)% of 
#Arnica associated pollinators are expected.


#model without Meligethes----
m_community2 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems), 
                       data = comb1, family = binomial)
summary(m_community2)

eff_community2 <- effect("log(Stems)",m_community2, xlevels = 50)  
eff.plot(eff_community2, plotdata = T,
         ylab = "Proportion of Arnica-associated pollinators",
         xlab = "Population size Arnica (Nr Stems)",
         main = "community composition without Meligethes",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
#qqnorm(resid(m_community)) #qqplot looks good
#qqnorm(resid(m_community2))

residuals_community2 <- simulateResiduals(fittedModel = m_community2)
plot(residuals_community2)
testOutliers(residuals_community2)
#residual vs predicted significant deviations


###effect sizes2----
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