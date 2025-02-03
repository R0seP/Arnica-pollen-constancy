#This script selects a model/models to examines whether the different pollinator 
#species use Arnica to different amounts.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
comb_all2 <- read.csv("comb_all2.csv", h = T)
comb_all2$Species <- as.factor(comb_all2$Species)

#define special_log function that sets the log of 0 to 0 for transforming the 
  #Nr of Arnica pollen with natural log later
special_log <- function(x) {
  sapply(x, function(element) {
    if (is.na(element)) {
      return(NA)
    } else if (element == 0) {
      return(0)
    } else {
      return(log(element))
    }
  })
}

#data visualization----
ggplot(comb_all2, aes(fill=Group, y=P_ASTE.Arnica.montana, x=Species)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(legend.text = element_text(size = 12), # increase legend text size
        legend.title = element_text(size = 14), # increase legend title size
        axis.text = element_text(size = 12), # increase axis text size
        axis.title = element_text(size = 14)) + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Species", y = "% Arnica montana pollen carried")

#relevel species----
#mean & median numbers of log of Arnica pollen per species and then the overall mean
comb_all2$log_NrArnica <- special_log(comb_all2$Nr_Arnica) #create column with log of Nr of Arncia where log(0) = 0

species_means <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(log_NrArnica, na.rm = TRUE))

#mean
overall_mean <- mean(na.omit(comb_all2$log_NrArnica))
overall_mean # = 3.18381

closest_value <- species_means %>%
  filter(abs(mean_value - overall_mean) == min(abs(mean_value - overall_mean)))

print(closest_value)
#Andrena sp carries on average a log number of pollen closest to overall average,
#use Andrena sp as baseline? (default anyways)

#median
overall_median <- median(na.omit(comb_all2$log_NrArnica))
overall_median # = 2.8

closest_value2 <- species_means %>%
  filter(abs(mean_value - overall_median) == min(abs(mean_value - overall_median)))

print(closest_value2)
#Bombus terrestris closest to median of proportion Arnica carried,
#use Bombus terrestris as baseline?


#model selection----
m1 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ (Stems + Group + Species)^2
               + (1|Site), family = binomial, 
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))) 
#control to increase the number of iterations to allow model to converge!
summary(m1)
check_overdispersion(m1)
#no overdispersion but rank deficient, multicollinearity
check_collinearity(m1)
#high correlation detected, maybe not ideal model

m2 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Species + Stems * Group
              + (1|Site), family = binomial,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m2)
check_overdispersion(m2)
check_collinearity(m2)
#no longer rank deficient. All interactions significant.
#still multicollinearity issues

m2.2 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Species + Stems * Group
                + (1|Site), family = binomial,
                data = comb_all2,
                control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m2.2)
check_overdispersion(m2.2)
check_collinearity(m2.2)
#no multicollinearity issues anymore

m3 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Species + log(Stems) * Group
              + (1|Site), family = binomial,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m3)
check_overdispersion(m3)
check_collinearity(m3)
#multicollinearity issues

mlist = list(m1, m2, m2.2, m3)
AICTab = AIC(m1, m2, m2.2, m3) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#for binomial models: m1 < m2 < m2.2 < m3, 
#all but m2.2 multicollinearity issues

#try poisson model because fit issues with binomial model detected
m4 <- glmmTMB(Nr_Arnica ~ (Species + Stems + Group)^2 + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m4)

check_overdispersion(m4)
#overdispersed
check_collinearity(m4)
#species and group collinear

m4.2 <- glmmTMB(Nr_Arnica ~ (Species + Stems + Group)^2 + offset(log(nPoll))
                + (1|Site), family = nbinom2,
                data = comb_all2,
                control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

summary(m4.2)

check_overdispersion(m4.2)
#underdispersed
check_collinearity(m4.2)
#species and group highly collinear, stems moderate correlation

m4.3 <- glmmTMB(Nr_Arnica ~ (Species + Stems + Group)^2 + offset(log(nPoll))
                + (1|Site), family = nbinom1,
                data = comb_all2,
                control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

summary(m4.3)

check_overdispersion(m4.3)
#no overdispersion
check_collinearity(m4.3)
#species high correlation, group moderate correlation

m5 <- glmmTMB(Nr_Arnica ~ log(Stems) * Species + log(Stems) * Group + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m5)

check_overdispersion(m5)
#overdispersed
check_collinearity(m5)
#species and stems, group and stems collinear

m6 <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m6)
check_overdispersion(m6)
#overdispersed
check_collinearity(m6)
#no collinearity

#negative binomial models because of overdispersion
m7 <- glmmTMB(Nr_Arnica ~ Stems * Species + Stems * Group + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m7)
check_overdispersion(m7)
#no overdipsersion anymore
check_collinearity(m7)
#species and stems collinear

m8 <- glmmTMB(Nr_Arnica ~ log(Stems) * Species + log(Stems) * Group + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m8)
check_overdispersion(m8)
#no overdispersion
check_collinearity(m8)
#species and group, stems and group, stems and species collinear

m9 <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m9)
check_overdispersion(m9)
#underdispersed
check_collinearity(m9)
#no multicollinearity

#try zero-inflated model to handle underdispersion
m9.2 <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                data = comb_all2, family = nbinom2, ziformula = ~1,
                control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m9.2)
check_overdispersion(m9.2)
#still underdispersed
check_collinearity(m9.2)
#no collinearity

#generalized poisson also still underdispersed, conway-maxwell distribution
#does not converge

#try to model variance as a linear function (nbinom1) instead to handle underdispersion
m9.3 <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll))
              + (1|Site), family = nbinom1,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m9.3)
check_overdispersion(m9.3)
#no overdispersion!
check_collinearity(m9.3)
#no collinearity!

#try with nPoll as predictor, not as offset
m10 <- glmmTMB(Nr_Arnica ~ (Stems * Species + Stems * Group) * log(nPoll)
               + (1|Site), family = nbinom2,
               data = comb_all2,
               control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m10)
check_overdispersion(m10)
#no overdispersion
check_collinearity(m10)
#diverse multicollinearity, among others with n(Poll)

m11 <- glmmTMB(Nr_Arnica ~ (Species + Stems * Group) * log(nPoll)
               + (1|Site), family = nbinom2,
               data = comb_all2,
               control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m11)
check_overdispersion(m11)
#underdispersed
check_collinearity(m11)
#still diverse multicollinearity, especially with n(Poll)


mlist = list(m4, m4.2, m4.3, m5, m6, m7, m8, m9, m9.2, m9.3, m10, m11)
AICTab = AIC(m4, m4.2, m4.3, m5, m6, m7, m8, m9, m9.2, m9.3, m10, m11) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#ranked m8 < m4.2 < m9 < m9.2 < m11 < m7 < m10 < m4.3 < m9.3 < m5 < m4 < m6.  
#All but m9.3 are over- or underdispersed or have multicollinearity issues
#m4.3 has some collinearity but is also not over-or underdisersed

r.squaredGLMM(m1) #does not currently work
r.squaredGLMM(m2)
r.squaredGLMM(m2.2)
r.squaredGLMM(m3)
r.squaredGLMM(m4) #does not currently work
r.squaredGLMM(m4.2) #does not currently work
r.squaredGLMM(m4.3) #does not currently work
r.squaredGLMM(m5)
r.squaredGLMM(m6)
r.squaredGLMM(m7)
r.squaredGLMM(m8)
r.squaredGLMM(m9)
r.squaredGLMM(m9.2)
r.squaredGLMM(m9.3)
r.squaredGLMM(m10)
r.squaredGLMM(m11)
#negative binomial models much lower r2.
#m9.3 (n.b. which does not have issues) 3rd lowest r2 (after other "9" models)
#m2.2 (b. which does not have issues) lowest r2 of binomial models, still very high

#model binomial----
m_species_binomial <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Species + Stems * Group
                     + (1|Site), family = binomial,
                     data = comb_all2,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_binomial)

eff_species1 <- effect("Stems",m_species_binomial, xlevels = 50)  
eff.plot(eff_species1, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "binomial model",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
#qqnorm(resid(m_species))
hist(resid(m_species_binomial)) #residual normality seems fine

residuals_binomial <- simulateResiduals(fittedModel = m_species_binomial)
plot(residuals_binomial)
testOutliers(residuals_binomial)
#all significant, model fit issues?

#model negative binomial----
#(m9.3)
m_species_negbi <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                           data = comb_all2, family = nbinom1,
                           control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_negbi)

eff_species2 <- effect(c("Stems"),m_species_negbi, xlevels = 50)
eff.plot(eff_species2, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "negative binomial model",
         ylim.data = T, overlay = F, 
         col.data = 3)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_species_negbi) #no overdipersion
check_collinearity(m_species_negbi) #no highly correlated predictors

hist(resid(m_species_negbi)) #distribution of residuals has tail to negative values

residuals_negbi <- simulateResiduals(fittedModel = m_species_negbi)
plot(residuals_negbi)
testOutliers(residuals_negbi)


#try with model with all interactions but with some multicollinearity (m4.3):
m_species_negbi2 <- glmmTMB(Nr_Arnica ~ (Species + Stems + Group)^2 + offset(log(nPoll))
                            + (1|Site), family = nbinom1,
                            data = comb_all2,
                            control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_negbi2)

eff_species3 <- effect(c("Stems"),m_species_negbi2, xlevels = 50)
eff.plot(eff_species3, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "negative binomial model",
         ylim.data = T, overlay = F, 
         col.data = 3)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_species_negbi2) #no overdipersion
check_collinearity(m_species_negbi2) #species and group high VIF, stems moderately high vif

hist(resid(m_species_negbi2)) #distribution of residuals has tail to negative values

residuals_negbi2 <- simulateResiduals(fittedModel = m_species_negbi2)
plot(residuals_negbi2)
testOutliers(residuals_negbi2)
#outlier test non-significant, residuals vs predicted looks strange, KS test 
#significant. Pred. vs resid. looks approximately similarly good to nb model
#1 (without interactions), chose model without collinearity (9.3, first negative
#binomial model)
