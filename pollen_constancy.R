#this script looks into flower constancy in my samples of pollinators on Arnica.

setwd("C:/Users/sohe1/documents/Master General Biology/Master_Thesis/R")
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
comb_all2$Nr_Arnica <- round(comb_all2$Nr_Arnica)
comb_all2$Nr_Not.Arnica <- round(comb_all2$Nr_Not.Arnica)
comb_flower <- subset(comb_all2, Group == "flower")   #n = 240
comb_area <- subset(comb_all2, Group == "area")     #n = 209


#general model----
#binomial model
m_pollen_gen <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                          (1|Site) + (1|Species), family = binomial, 
                        na.action = na.omit, data = comb_all2)
summary(m_pollen_gen)
#significant interactions between Stems and Group, therefore split the data set
#from here on in pollinators caught on Arnica (flower) and in the area (area)

check_overdispersion(m_pollen_gen)
#after rounding Nr Arnica no overdispersion anymore


#negative binomial model with offset
m_pollen_gen2 <- glmmTMB(Nr_Arnica ~ Stems * Group +
                          (1|Site) + (1|Species) + offset(log(nPoll)), 
                        family = nbinom2, 
                        na.action = na.omit, 
                        data = comb_all2)
summary(m_pollen_gen2)
check_overdispersion(m_pollen_gen2)
#now underdispersion

par(mfrow = c(1,2))
eff_pollen_gen1 <- effect("Stems",m_pollen_gen, xlevels = 50)  
eff.plot(eff_pollen_gen1, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "all pollinators, binomial model",
         ylim.data = T, overlay = F, col.data = 3) 

eff_pollen_gen2 <- effect("Stems",m_pollen_gen2, xlevels = 50)  
eff.plot(eff_pollen_gen2, plotdata = T,
         ylab = "Number of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "all pollinators, negative binomial model",
         ylim.data = T, overlay = F, col.data = 3)

#try model with weights instead of offset (negative binomial model)
m_pollen_gen2.2 <- glmmTMB(Nr_Arnica ~ Stems * Group +
                             (1|Site) + (1|Species), 
                           weights = nPoll,
                           family = nbinom2, 
                           na.action = na.omit, 
                           data = comb_all2)
summary(m_pollen_gen2.2)
check_overdispersion(m_pollen_gen2.2)
#also now underdispersed

eff_pollen_gen2 <- effect("Stems",m_pollen_gen2, xlevels = 50)  
eff.plot(eff_pollen_gen2, plotdata = T,
         ylab = "Number of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "n.b. model with offset",
         ylim.data = T, overlay = F, col.data = 3)
eff_pollen_gen2.2 <- effect("Stems",m_pollen_gen2.2, xlevels = 50)  
eff.plot(eff_pollen_gen2.2, plotdata = T,
         ylab = "Number of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "n.b. model with weights",
         ylim.data = T, overlay = F, col.data = 3)

#see whether models with log(Stems) better
#binomial model with log(Stems)
m_pollen_gen3 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group +
                          (1|Site) + (1|Species), family = binomial, 
                        na.action = na.omit, data = comb_all2)
summary(m_pollen_gen3)
check_overdispersion(m_pollen_gen3)
#after rounding Nr Arnica no overdispersion anymore

eff_pollen_gen3 <- effect("log(Stems)",m_pollen_gen3, xlevels = 50)  
eff.plot(eff_pollen_gen3, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "all pollinators, b. model",
         ylim.data = T, overlay = F, col.data = 3)

#negative binomial model with offset and log(Stems)
m_pollen_gen4 <- glmmTMB(Nr_Arnica ~ log(Stems) * Group +
                           (1|Site) + (1|Species) + offset(log(nPoll)), 
                         family = nbinom2,
                         na.action = na.omit, 
                         data = comb_all2)
summary(m_pollen_gen4)
check_overdispersion(m_pollen_gen4)
#also underdispersed

eff_pollen_gen4 <- effect("log(Stems)",m_pollen_gen4, xlevels = 50)  
eff.plot(eff_pollen_gen4, plotdata = T,
         ylab = "Number of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "all pollinators, n.b. model",
         ylim.data = T, overlay = F, col.data = 3)

par(mfrow=c(1,1))

#model with nPoll as predictor (not weights or offset)
m_pollen_gen5 <- glmmTMB(Nr_Arnica ~ (Stems + Group + log(nPoll))^2 +
                           (1|Site) + (1|Species) , 
                         family = nbinom2,
                         na.action = na.omit, 
                         data = comb_all2)
summary(m_pollen_gen5)
check_overdispersion(m_pollen_gen5)
#underdispersed

m_pollen_gen5b <- glmmTMB(Nr_Arnica ~ Stems * Group + Group * log(nPoll) +
                            (1|Site) + (1|Species) , 
                          family = nbinom2,
                          na.action = na.omit, 
                          data = comb_all2)
summary(m_pollen_gen5b)

m_pollen_gen5c <- glmmTMB(Nr_Arnica ~ Stems + Group * log(nPoll) +
                            (1|Site) + (1|Species) , 
                          family = nbinom2,
                          na.action = na.omit, 
                          data = comb_all2)
summary(m_pollen_gen5c)
check_overdispersion(m_pollen_gen5c)
#also underdisersed



eff_pollen_gen5 <- effect("Stems",m_pollen_gen5, xlevels = 50)  
eff.plot(eff_pollen_gen5, plotdata = T,
         ylab = "Number of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "all pollinators, n.b. model with nPoll as predictor",
         ylim.data = T, overlay = F, col.data = 3)





mlist = list(m_pollen_gen, m_pollen_gen2, m_pollen_gen2.2, m_pollen_gen3, 
             m_pollen_gen4, m_pollen_gen5, m_pollen_gen5b, m_pollen_gen5c)
AICTab = AIC(m_pollen_gen, m_pollen_gen2, m_pollen_gen2.2, m_pollen_gen3, 
             m_pollen_gen4, m_pollen_gen5, m_pollen_gen5b, m_pollen_gen5c) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#negative binomial models with offset ranked higher than binomial ranked higher 
#than negative binomial model with weights



#model selection----
###flower----
m1 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                (1|Site) + (1|Species), family = binomial, 
              na.action = na.omit, data = comb_flower)
m2 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) +
                (1|Site) + (1|Species), family = binomial, 
              na.action = na.omit, data = comb_flower)

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
#model m2 with log also higher r2, use that model!

###area----
m1 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                (1|Site) + (1|Species), family = binomial, 
              na.action = na.omit, data = comb_area)
m2 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) +
                (1|Site) + (1|Species), family = binomial, 
              na.action = na.omit, data = comb_area)

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
#model m2 with log also higher r2, use that model!

#models----
#modelling the impact of Arnica population size on Arnica pollen carried by 
#pollinators caught on flowers and in area separately

###flower model----
#for pollinators caught on flowers 
m_pollen_flower <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) +
                             (1|Site) + (1|Species), family = binomial, 
                           na.action = na.omit, data = comb_flower)

summary(m_pollen_flower)
check_overdispersion(m_pollen_flower)
#no overdispersion

eff_constancy_flower <- effect("log(Stems)",m_pollen_flower, xlevels = 50)  
eff.plot(eff_constancy_flower, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "pollinators caught on flower",
         ylim.data = T, overlay = F, col.data = 3) 

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_pollen_flower))  
hist(resid(m_pollen_flower))
#looks very good

resid_jittered <- resid(m_pollen_flower) + rnorm(length(resid(m_pollen_flower)), mean = 0, sd = 1e-6)
ks.test(resid_jittered, "pnorm", mean = mean(resid_jittered), sd = sd(resid_jittered))
#no significance, so no deviation of the residuals from a normal distribution

residuals_pollen_flower <- simulateResiduals(fittedModel = m_pollen_flower, re.form = ~0)
plot(residuals_pollen_flower)
testOutliers(residuals_pollen_flower)
#KS and outlier test significant, residual vs predicted significant

###area model----
###for pollinators caught in the area
m_pollen_area <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) +
                           (1|Site) + (1|Species), family = binomial, 
                         na.action = na.omit, data = comb_area)

summary(m_pollen_area)

eff_constancy_area <- effect("log(Stems)",m_pollen_area, xlevels = 50)
eff.plot(eff_constancy_area, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "pollinators caught in area",
         ylim.data = T, overlay = F, col.data = 3)

m_pollen_area <- glmmTMB(Nr_Arnica ~ Stems * log(nPoll) +
                           (1|Site) + (1|Species) , 
                         family = nbinom2,
                         na.action = na.omit, 
                         data = comb_area)
summary(m_pollen_area)
check_overdispersion(m_pollen_area)

m_pollen_area2 <- glmmTMB(Nr_Arnica ~ Stems + log(nPoll) +
                            (1|Site) + (1|Species) , 
                          family = nbinom2,
                          na.action = na.omit, 
                          data = comb_area)
summary(m_pollen_area2)
check_overdispersion(m_pollen_area2)


#test if model assumptions are met and test model for fit:
qqnorm(resid(m_pollen_area))  #qq looks iffy
hist(resid(m_pollen_area))  #only very small residuals, definitely no normality

ks.test(resid(m_pollen_area), "pnorm", mean = mean(resid(m_pollen_area)), sd = sd(resid(m_pollen_area)))
#significant, deviation of the residuals from a normal distribution, model may 
#not be a good fit?

residuals_pollen_area <- simulateResiduals(fittedModel = m_pollen_area)
plot(residuals_pollen_area)
testOutliers(residuals_pollen_area)
#does not work at the moment, produces infinite values, try again!

#effect sizes----
####flower----
#r-squared value:
r.squaredGLMM(m_pollen_flower)

#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_flower <- data.frame(Stems = c(10, 100, 500))
predictions_flower <- predict(m_pollen_flower, newdata = pred_data_flower, type = "response", re.form = NA, se.fit = T)

pred_data_flower$Pred.Arnica <- predictions_flower$fit
pred_data_flower$Pred.Arnica_SE <- predictions_flower$se.fit
pred_data_flower
#for 10 stems, 32.0 (+/- 7.0)% of Arnica pollen are expected, for 100 stems
#44.1 (+/- 5.7)% of Arnica pollen, and for 500 stems 53.1 (+/- 8.1)% of 
#Arnica pollen are expected to be carried by a pollinator caught on Arnica, 
#averaged over all sites.

####area----
r.squaredGLMM(m_pollen_area)

#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_area <- data.frame(Stems = c(10, 100, 500))
predictions_area <- predict(m_pollen_area, newdata = pred_data_area, type = "response", re.form = NA, se.fit = T)

pred_data_area$Pred.Arnica <- predictions_area$fit
pred_data_area$Pred.Arnica_SE <- predictions_area$se.fit
pred_data_area
#for 10 stems, 1.6 (+/- 2.1)% of Arnica pollen are expected, for 100 stems
#0.5 (+/- 0.5)% of Arnica pollen, and for 500 stems 0.2 (+/- 0.3)% of 
#Arnica pollen are expected to be carried by a pollinator caught in the area, 
#averaged over all sites.


#models subgroups----
#modelling the impact of Arnica population size on Arnica pollen in subgroups of
#pollinators associated with Arnica and pollinators associated with the area

#this part of the analysis is not currently used!
#unsure about justification to split data set and about expectations for each subgroup!

#all Arnica associated species
m_Arnica_associated <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                                 (1|Site), family = binomial, 
                               na.action = na.omit, data = Arnica_polli)
summary(m_Arnica_associated)
plot(P_ASTE.Arnica.montana ~ Stems, data = Arnica_polli, main = "Arnica associated pollinators")
plot(effect("Stems",m_Arnica_associated),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Arnica associated pollinators")


#all area associated species
m_Area_associated <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                               (1|Site), family = binomial, 
                             na.action = na.omit, data = Area_polli)
summary(m_Area_associated)
plot(P_ASTE.Arnica.montana ~ Stems, data = Area_polli, main = "Area associated pollinators")
plot(effect("Stems",m_Area_associated),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Area associated pollinators")


#Arnica associated species caught on flower
subset_Arnica_polli <- subset(Arnica_polli, Group == "flower")
m_Arnica_associated_flower <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                        (1|Site), family = binomial, 
                                      na.action = na.omit, data = subset_Arnica_polli)
summary(m_Arnica_associated_flower)

plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Arnica_polli, main = "Arnica associated pollinators flower")
plot(effect("Stems",m_Arnica_associated_flower),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Arnica associated pollinators flower")

#Arnica associated species caught in the area
subset_Arnica_polli2 <- subset(Arnica_polli, Group == "area")
m_Arnica_associated_area <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                      (1|Site), family = binomial, 
                                    na.action = na.omit, data = subset_Arnica_polli2)
summary(m_Arnica_associated_area)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Arnica_polli2, main = "Arnica associated pollinators area")
plot(effect("Stems",m_Arnica_associated_area),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Arnica associated pollinators area")


#Area associated pollinators caught on flower
subset_Area_polli <- subset(Area_polli, Group == "flower")
m_Area_associated_flower <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                      (1|Site), family = binomial, 
                                    na.action = na.omit, data = subset_Area_polli)
summary(m_Area_associated_flower)

plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Area_polli, main = "Area associated pollinators flower")
plot(effect("Stems",m_Area_associated_flower),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Area associated pollinators flower")

#Area associated pollinators caught in area
subset_Area_polli2 <- subset(Area_polli, Group == "area")
m_Area_associated_area <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                    (1|Site), family = binomial, 
                                  na.action = na.omit, data = subset_Area_polli2)
summary(m_Area_associated_area)

plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Area_polli2, main = "Area associated pollinators area")
plot(effect("Stems",m_Area_associated_area),
     las=2,
     ylab = "% Arnica carried",
     xlab = "Population size Arnica (Nr Stems)",
     main = "Area associated pollinators area")
