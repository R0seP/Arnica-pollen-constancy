#This script examines whether the different pollinator species use Arnica to 
#different amounts.

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

comb_all2$Species <- as.factor(comb_all2$Species)

#comb_imp_species <- comb_all2 %>%
 # semi_join(imp_species, by = c("Species" = "Species"))

#define special_log function that sets the log of 0 to 0 for transforming the 
  #Nr of Arnica pollen with natural log later
special_log <- function(x) {
  sapply(x, function(element) {
    if (is.na(element)) {
      return(NA)
    } else if (element == 0) {
      return(0)
    } else {
      return(round(log(element)))
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

#mean & median levels of Arnica per species----
# Calculate means for each species and then the overall mean
species_means <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE))

#mean
overall_mean <- mean(species_means$mean_value)
overall_mean # = 0.2477039

closest_value <- species_means %>%
  filter(abs(mean_value - overall_mean) == min(abs(mean_value - overall_mean)))

print(closest_value)
#Apis mellifera closest to mean of amount Arnica carried,
#use Apis mellifera as baseline?

#median
overall_median <- median(species_means$mean_value)
overall_median # = 0.1813998

closest_value2 <- species_means %>%
  filter(abs(mean_value - overall_median) == min(abs(mean_value - overall_median)))

print(closest_value2)
#Maniola jurtina sp closest to median of proportion Arnica carried,
#use Maniola jurtina as baseline?

#model selection----
comb_all2$Species <- relevel(comb_all2$Species, ref = "Apis mellifera")
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
#species and stems collinear

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


mlist = list(m4, m5, m6, m7, m8, m9, m9.2, m9.3, m10, m11)
AICTab = AIC(m4, m5, m6, m7, m8, m9, m9.2, m9.3, m10, m11) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#ranked m8 < m9 < m9.2 < m11 < m7 < m10 < m9.3 < m5 < m4 < m6.  
#All but m9.3 are over- or underdispersed or have multicollinearity issues

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m2.2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
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
m_species_negbi <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                           data = comb_all2, family = nbinom1,
                           control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_negbi)
check_overdispersion(m_species_negbi)

eff_species2 <- effect(c("Stems"),m_species_negbi, xlevels = 50)
eff.plot(eff_species2, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "negative binomial model",
         ylim.data = T, overlay = F, 
         col.data = 3)

#test if model assumptions are met and test model for fit:
#qqnorm(resid(m_species2))
hist(resid(m_species_negbi)) #distribution of residuals has tail to negative values

residuals_negbi <- simulateResiduals(fittedModel = m_species_negbi)
plot(residuals_negbi)
testOutliers(residuals_negbi)


#try with log of response:
m_species_negbi2 <- glmmTMB(special_log(Nr_Arnica) ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                            data = comb_all2, family = nbinom1,
                            control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_negbi2)

eff_species3 <- effect(c("Stems"),m_species_negbi2, xlevels = 50)
eff.plot(eff_species3, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "negative binomial model",
         ylim.data = T, overlay = F, 
         col.data = 3)

#test if model assumptions are met and test model for fit:
#qqnorm(resid(m_species2))
hist(resid(m_species_negbi2)) #distribution of residuals has tail to negative values

residuals_negbi2 <- simulateResiduals(fittedModel = m_species_negbi2)
plot(residuals_negbi2)
testOutliers(residuals_negbi2)
#outlier test non-significant, residuals vs predicted looks strange, KS test 
#significant. Seems worse than model that uses normal response!

#effect sizes----

#r-squared value:
r.squaredGLMM(m_species)

#predict the percentage of Arnica pollen for species in different groups for 10, 100, 500 stems
pred_data_species <- data.frame(Stems = c(rep(10,44),rep(100,44),rep(500,44)), 
                                Species = rep(imp_species$Species,6), 
                                Group = c(rep("area",22), rep("flower",22)))
predictions_species <- predict(m_species, newdata = pred_data_species, type = "response", re.form = NA, se.fit = T)

pred_data_species$Pred.Arnica <- predictions_species$fit
pred_data_species$Pred.Arnica_SE <- predictions_species$se.fit
pred_data_species

#models by species----

###Andrena sp----
Andrena <- subset(comb_all2, Species == "Andrena sp")

mAndrena <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ (log(Stems) * Group)
                    + (1|Site), family = binomial,
                    data = Andrena,
                    control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mAndrena)
summary(mAndrena)

eff_Andrena <- effect(c("log(Stems)"),mAndrena, xlevels = 50)
eff.plot(eff_Andrena, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Adrena sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Andrena <- simulateResiduals(fittedModel = mAndrena)
plot(residuals_Andrena)
testOutliers(residuals_Andrena)


###Apis mellifera----
Apis <- subset(comb_all2, Species == "Apis mellifera")

#binomial model
mApis <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ (log(Stems) * Group)
                    + (1|Site), family = binomial,
                    data = Apis,
                    control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mApis)
summary(mApis)

eff_Apis <- effect(c("log(Stems)"),mApis, xlevels = 50)
eff.plot(eff_Apis, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Apis mellifera",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Apis <- simulateResiduals(fittedModel = mApis)
plot(residuals_Apis)
testOutliers(residuals_Apis)

#negative binomial model
mApis2 <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                  + (1|Site), family = nbinom1,
                  data = Apis)

check_overdispersion(mApis2)
#underdispersion!
summary(mApis2)

eff_Apis2 <- effect(c("Stems"),mApis2, xlevels = 50)
eff.plot(eff_Apis2, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Apis mellifera",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Apis2 <- simulateResiduals(fittedModel = mApis2)
plot(residuals_Apis2)
testOutliers(residuals_Apis2)
#significant KS and significant res vs pred, binomial model much better fit here

###Bombus pascuorum----
Bombus_pascuorum <- subset(comb_all2, Species == "Bombus pascuorum")

mBombus_pascuorum <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems)
                 + (1|Site), family = binomial,
                 data = Bombus_pascuorum,
                 control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mBombus_pascuorum)
summary(mBombus_pascuorum)

eff_Bombus_pascuorum <- effect(c("log(Stems)"),mBombus_pascuorum, xlevels = 50)
eff.plot(eff_Bombus_pascuorum, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Bombus pascuorum",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Bombus_pascuorum <- simulateResiduals(fittedModel = mBombus_pascuorum)
plot(residuals_Bombus_pascuorum)
testOutliers(residuals_Bombus_pascuorum)


###Bombus ruderarius----
Bombus_ruderarius <- subset(comb_all2, Species == "Bombus ruderarius")

mBombus_ruderarius <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) 
                             + (1|Site), family = binomial,
                             data = Bombus_ruderarius,
                             control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mBombus_ruderarius)
summary(mBombus_ruderarius)

eff_Bombus_ruderarius <- effect(c("log(Stems)"),mBombus_ruderarius, xlevels = 50)
eff.plot(eff_Bombus_ruderarius, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Bombus ruderarius",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Bombus_ruderarius <- simulateResiduals(fittedModel = mBombus_ruderarius)
plot(residuals_Bombus_ruderarius)
testOutliers(residuals_Bombus_ruderarius)


###Bombus terrestris----
Bombus_terrestris <- subset(comb_all2, Species == "Bombus terrestris")

mBombus_terrestris <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group 
                              + (1|Site), family = binomial,
                              data = Bombus_terrestris,
                              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mBombus_terrestris)
summary(mBombus_terrestris)

eff_Bombus_terrestris <- effect(c("log(Stems)"),mBombus_terrestris, xlevels = 50)
eff.plot(eff_Bombus_terrestris, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Bombus terrestris",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Bombus_terrestris <- simulateResiduals(fittedModel = mBombus_terrestris)
plot(residuals_Bombus_terrestris)
testOutliers(residuals_Bombus_terrestris)


###Coenonympha pamphilus----
Coenonympha <- subset(comb_all2, Species == "Coenonympha pamphilus")

mCoenonympha <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems)
                              + (1|Site), family = binomial,
                              data = Coenonympha,
                              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mCoenonympha)
summary(mCoenonympha)

eff_Coenonympha <- effect(c("log(Stems)"),mCoenonympha, xlevels = 50)
eff.plot(eff_Coenonympha, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Coenonympha pamphilus",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Coenonympha <- simulateResiduals(fittedModel = mCoenonympha)
plot(residuals_Coenonympha)
testOutliers(residuals_Coenonympha)


###Dasytes niger----
Dasytes <- subset(comb_all2, Species == "Dasytes niger")

mDasytes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) + Group
                        + (1|Site), family = binomial,
                        data = Dasytes,
                        control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mDasytes)
summary(mDasytes)

eff_Dasytes <- effect(c("log(Stems)"),mDasytes, xlevels = 50)
eff.plot(eff_Dasytes, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Dasytes niger",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Dasytes <- simulateResiduals(fittedModel = mDasytes)
plot(residuals_Dasytes)
testOutliers(residuals_Dasytes)


###Empis livida----
Empis_livida <- subset(comb_all2, Species == "Empis livida")

mEmpis_livida <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) 
                    + (1|Site), family = binomial,
                    data = Empis_livida,
                    control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mEmpis_livida)
summary(mEmpis_livida)

eff_Empis_livida <- effect(c("log(Stems)"),mEmpis_livida, xlevels = 50)
eff.plot(eff_Empis_livida, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Empis livida",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Empis_livida <- simulateResiduals(fittedModel = mEmpis_livida)
plot(residuals_Empis_livida)
testOutliers(residuals_Empis_livida)


###Empis tessellata----
Empis_tessellata <- subset(comb_all2, Species == "Empis tessellata")

mEmpis_tessellata <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group
                         + (1|Site), family = binomial,
                         data = Empis_tessellata,
                         control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))


check_overdispersion(mEmpis_tessellata)
summary(mEmpis_tessellata)

eff_Empis_tessellata <- effect(c("log(Stems)"),mEmpis_tessellata, xlevels = 50)
eff.plot(eff_Empis_tessellata, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Empis tessellata",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Empis_tessellata <- simulateResiduals(fittedModel = mEmpis_tessellata)
plot(residuals_Empis_tessellata)
testOutliers(residuals_Empis_tessellata)


###Eristalis sp----
Eristalis <- subset(comb_all2, Species == "Eristalis sp")

Eristalis <- Eristalis[Eristalis$Nr_Arnica != 0, ]

mEristalis <- glmmTMB(round(log(Nr_Arnica)) ~ Stems * Group + offset(log(nPoll))
                              + (1|Site), family = nbinom2,
                             data = Eristalis)

check_overdispersion(mEristalis)
#no overdispersion with nbinom2
summary(mEristalis)

eff_Eristalis <- effect(c("Stems"),mEristalis, xlevels = 50)
eff.plot(eff_Eristalis, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Eristalis sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Eristalis <- simulateResiduals(fittedModel = mEristalis)
plot(residuals_Eristalis)
testOutliers(residuals_Eristalis)


###Eupeodes corollae----
Eupeodes <- subset(comb_all2, Species == "Eupeodes corollae")

Eupeodes_no0 <- Eupeodes[Eupeodes$Nr_Arnica != 0, ]

#negative binomial model
mEupeodes1 <- glmmTMB(round(log(Nr_Arnica)) ~ Stems * Group + offset(log(nPoll))
                     + (1|Site), family = nbinom2,
                     data = Eupeodes_no0)

check_overdispersion(mEupeodes1)
summary(mEupeodes1)

eff_Eupeodes1 <- effect(c("Stems"),mEupeodes1, xlevels = 50)
eff.plot(eff_Eupeodes1, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Eupeodes corollae",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Eupeodes1 <- simulateResiduals(fittedModel = mEupeodes1)
plot(residuals_Eupeodes1)
testOutliers(residuals_Eupeodes1)
#residuals vs predicted significant

#negative binomial model without log on reponse:
mEupeodes1.2 <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                      + (1|Site), family = nbinom2,
                      data = Eupeodes)

check_overdispersion(mEupeodes1.2)
summary(mEupeodes1.2)

eff_Eupeodes1.2 <- effect(c("Stems"),mEupeodes1.2, xlevels = 50)
eff.plot(eff_Eupeodes1.2, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Eupeodes corollae",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Eupeodes1.2 <- simulateResiduals(fittedModel = mEupeodes1.2)
plot(residuals_Eupeodes1.2)
testOutliers(residuals_Eupeodes1.2)
#residuals vs predicted significant, maybe slightly worse than for model with 
  #log response, but uses all data points

#binomial model
mEupeodes2 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group
                      + (1|Site), family = binomial,
                      data = Eupeodes)

check_overdispersion(mEupeodes2)
summary(mEupeodes2)

eff_Eupeodes2 <- effect(c("Stems"),mEupeodes2, xlevels = 50)
eff.plot(eff_Eupeodes2, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Eupeodes corollae",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Eupeodes2 <- simulateResiduals(fittedModel = mEupeodes2)
plot(residuals_Eupeodes2)
testOutliers(residuals_Eupeodes2)
#KS test and residuals vs predicted significant
#res vs pred same weird directional tendency as Eristalis and Eupeodes negative
  #binomial model, but worse than for Eupeodes nb model

###Helophilus pendulus----
Helophilus <- subset(comb_all2, Species == "Helophilus pendulus")

mHelophilus <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group
                     + (1|Site), family = binomial,
                     data = Helophilus,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mHelophilus)
summary(mHelophilus)

eff_Helophilus <- effect(c("log(Stems)"),mHelophilus, xlevels = 50)
eff.plot(eff_Helophilus, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Helophilus pendulus",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Helophilus <- simulateResiduals(fittedModel = mHelophilus)
plot(residuals_Helophilus)
testOutliers(residuals_Helophilus)


###Lasioglossum sp----
Lasioglossum <- subset(comb_all2, Species == "Lasioglossum sp")

mLasioglossum <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group
                       + (1|Site), family = binomial,
                       data = Lasioglossum,
                       control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mLasioglossum)
summary(mLasioglossum)

eff_Lasioglossum <- effect(c("log(Stems)"),mLasioglossum, xlevels = 50)
eff.plot(eff_Lasioglossum, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Lasioglossum sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Lasioglossum <- simulateResiduals(fittedModel = mLasioglossum)
plot(residuals_Lasioglossum)
testOutliers(residuals_Lasioglossum)


###Maniola jurtina----
Maniola <- subset(comb_all2, Species == "Maniola jurtina")

mManiola <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) + Group
                         + (1|Site), family = binomial,
                         data = Maniola,
                         control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mManiola)
summary(mManiola)

eff_Maniola <- effect(c("log(Stems)"),mManiola, xlevels = 50)
eff.plot(eff_Maniola, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Maniola jurtina",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Maniola <- simulateResiduals(fittedModel = mManiola)
plot(residuals_Maniola)
testOutliers(residuals_Maniola)


###Meligethes sp----
Meligethes <- subset(comb_all2, Species == "Meligethes sp")

mMeligethes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) 
                    + (1|Site), family = binomial,
                    data = Meligethes,
                    control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mMeligethes)
summary(mMeligethes)

eff_Meligethes <- effect(c("log(Stems)"),mMeligethes, xlevels = 50)
eff.plot(eff_Meligethes, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Meligethes sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Meligethes <- simulateResiduals(fittedModel = mMeligethes)
plot(residuals_Meligethes)
testOutliers(residuals_Meligethes)


###Merodon equestris----
Merodon <- subset(comb_all2, Species == "Merodon equestris")

mMerodon <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems)* Group
                       + (1|Site), family = binomial,
                       data = Merodon,
                       control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

check_overdispersion(mMerodon)
summary(mMerodon)

eff_Merodon <- effect(c("log(Stems)"),mMerodon, xlevels = 50)
eff.plot(eff_Merodon, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Merodon equestris",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Merodon <- simulateResiduals(fittedModel = mMerodon)
plot(residuals_Merodon)
testOutliers(residuals_Merodon)


###Nomada sp----
Nomada <- subset(comb_all2, Species == "Nomada sp")

mNomada <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems)* Group
                    + (1|Site), family = binomial,
                    data = Nomada)

check_overdispersion(mNomada)
summary(mNomada)

eff_Nomada <- effect(c("log(Stems)"),mNomada, xlevels = 50)
eff.plot(eff_Nomada, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Nomada sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Nomada <- simulateResiduals(fittedModel = mNomada)
plot(residuals_Nomada)
testOutliers(residuals_Nomada)


###Ochlodes sylvanus----
Ochlodes <- subset(comb_all2, Species == "Ochlodes sylvanus")

mOchlodes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems)
                   + (1|Site), family = binomial,
                   data = Ochlodes)

check_overdispersion(mOchlodes)
summary(mOchlodes)

eff_Ochlodes <- effect(c("log(Stems)"),mOchlodes, xlevels = 50)
eff.plot(eff_Ochlodes, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Ochlodes sylvanus",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Ochlodes <- simulateResiduals(fittedModel = mOchlodes)
plot(residuals_Ochlodes)
testOutliers(residuals_Ochlodes)


###Oedemera sp----
Oedemera <- subset(comb_all2, Species == "Oedemera sp")

mOedemera <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group
                     + (1|Site), family = binomial,
                     data = Oedemera)

check_overdispersion(mOedemera)
summary(mOedemera)

eff_Oedemera <- effect(c("log(Stems)"),mOedemera, xlevels = 50)
eff.plot(eff_Oedemera, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Oedemera sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Oedemera <- simulateResiduals(fittedModel = mOedemera)
plot(residuals_Oedemera)
testOutliers(residuals_Oedemera)


###Phyllopertha horticola----
Phyllopertha <- subset(comb_all2, Species == "Phyllopertha horticola")

mPhyllopertha <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) + Group
                     + (1|Site), family = binomial,
                     data = Phyllopertha)

check_overdispersion(mPhyllopertha)
summary(mPhyllopertha)

eff_Phyllopertha <- effect(c("log(Stems)"),mPhyllopertha, xlevels = 50)
eff.plot(eff_Phyllopertha, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Phyllopertha horticola",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Phyllopertha <- simulateResiduals(fittedModel = mPhyllopertha)
plot(residuals_Phyllopertha)
testOutliers(residuals_Phyllopertha)


###Sphaerophoria sp----
Sphaerophoria <- subset(comb_all2, Species == "Sphaerophoria sp")

mSphaerophoria <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) + Group
                         + (1|Site), family = binomial,
                         data = Sphaerophoria)

check_overdispersion(mSphaerophoria)
summary(mSphaerophoria)

eff_Sphaerophoria <- effect(c("log(Stems)"),mSphaerophoria, xlevels = 50)
eff.plot(eff_Sphaerophoria, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Sphaerophoria sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Sphaerophoria <- simulateResiduals(fittedModel = mSphaerophoria)
plot(residuals_Sphaerophoria)
testOutliers(residuals_Sphaerophoria)


###Stenurella melanura----
Stenurella <- subset(comb_all2, Species == "Stenurella melanura")

mStenurella <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Group
                          + (1|Site), family = binomial,
                          data = Stenurella)

check_overdispersion(mStenurella)
summary(mStenurella)

eff_Stenurella <- effect(c("log(Stems)"),mStenurella, xlevels = 50)
eff.plot(eff_Stenurella, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Stenurella mealnura",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Stenurella <- simulateResiduals(fittedModel = mStenurella)
plot(residuals_Stenurella)
testOutliers(residuals_Stenurella)
