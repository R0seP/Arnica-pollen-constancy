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
#use Andrena sp as baseline?

#median
overall_median <- median(na.omit(comb_all2$Nr_Arnica))
overall_median # = 16

closest_value2 <- species_means %>%
  filter(abs(mean_value - overall_median) == min(abs(mean_value - overall_median)))

print(closest_value2)
#Eueodes corollae sp closest to median of proportion Arnica carried,
#use Eupeodes corollae as baseline?

#model selection----
comb_all2$Species <- relevel(comb_all2$Species, ref = "Nomada sp")
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
check_overdispersion(m_species_negbi) #no overdipersion
check_collinearity(m_species_negbi) #no highly correlated predictors

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


#try with model with all interactions but with some multicollinearity (m4.3):
m_species_negbi2 <- glmmTMB(Nr_Arnica ~ (Species + Stems + Group)^2 + offset(log(nPoll))
                            + (1|Site), family = nbinom1,
                            data = comb_all2,
                            control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species_negbi2)
check_overdispersion(m_species_negbi2) #no overdipersion
check_collinearity(m_species_negbi2) #species and group high VIF, stems moderately high vif

eff_species3 <- effect(c("Stems"),m_species_negbi2, xlevels = 50)
eff.plot(eff_species3, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
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
#significant. Pred. vs resid. looks approximately similarly good to nb model
#1 (without interactions), chose model without collinearity

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


#relevel species----
m_species <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                     data = comb_all2, family = nbinom1,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

#extract the species with the average model coefficient:
coefficients <- as.data.frame(summary(m_species)$coefficients$cond) #extract coefficients 

imp_species$Parameter <- coefficients[1:22,1] #add parameter estimates to data frame

imp_species_new <- imp_species #create new data frame
imp_species_new$Intercept <- rep(0,22)

imp_species_new[1,"Intercept"] <- imp_species_new[1,"Parameter"] #set intercept for Andrena to be same as parameter estimate for intercept

for(i in 1:21){
  imp_species_new[1+i, "Intercept"] <- imp_species[1+i,"Parameter"] + imp_species[1,"Parameter"]
} #calculate intercept for species by adding parameter to intercept parameter estimate

mean_intercept <- mean(imp_species_new$Intercept) #calculate mean intercept
mean_intercept # = -2.695824

closest_index <- which.min(abs(imp_species_new$Intercept - mean_intercept)) #find row with intercept closest to mean intercept
closest_row <- imp_species_new[closest_index, ]
print(closest_row)
#Stenurella melanura species with mean intercept, use as baseline?


#models by species----
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
