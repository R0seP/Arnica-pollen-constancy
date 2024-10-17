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
comb_all2$Nr_Arnica <- round(comb_all2$Nr_Arnica)
comb_all2$Nr_Not.Arnica <- round(comb_all2$Nr_Not.Arnica)
comb_all2$Species <- as.factor(comb_all2$Species)

#comb_imp_species <- comb_all2 %>%
 # semi_join(imp_species, by = c("Species" = "Species"))

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
overall_mean # = 0.2332081

closest_value <- species_means %>%
  filter(abs(mean_value - overall_mean) == min(abs(mean_value - overall_mean)))

print(closest_value)
#Maniola jurtina closest to mean of amount Arnica carried,
#use Maniola jurtina as baseline?

#median
overall_median <- median(species_means$mean_value)
overall_median # = 0.17553538

closest_value2 <- species_means %>%
  filter(abs(mean_value - overall_median) == min(abs(mean_value - overall_median)))

print(closest_value2)
#Andrena sp closest to median of proportion Arnica carried,
#use Andrena sp as baseline?

#model selection----
comb_all2$Species <- relevel(comb_all2$Species, ref = "Maniola jurtina")
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
#no longer rank deficient. All interactions significant.

m3 <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Species + log(Stems) * Group
              + (1|Site), family = binomial,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m3)

mlist = list(m1, m2, m3)
AICTab = AIC(m1, m2, m3) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#for binomial models: m1 < m3 < m2, but m1 multicollinearity issues, so m3 best

#try poisson model because fit issues with binomial model detected
m4 <- glmmTMB(Nr_Arnica ~ Stems * Species + Stems * Group + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m4)
#here, interaction Stems and Group not significant, so take out
check_overdispersion(m4)
#overdispersed

m5 <- glmmTMB(Nr_Arnica ~ log(Stems) * Species + log(Stems) * Group + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m5)
#here, interaction Stems and Group significant
check_overdispersion(m5)
#overdispersed

m6 <- glmmTMB(Nr_Arnica ~ Stems * Species + Group + offset(log(nPoll))
              + (1|Site), family = poisson,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m6)
check_overdispersion(m6)
#overdispersed


#negative binomial models because of overdispersion
m7 <- glmmTMB(Nr_Arnica ~ Stems * Species + Stems * Group + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m7)
check_overdispersion(m7)
#no overdipsersion anymore, interaction stems and group non significant

m8 <- glmmTMB(Nr_Arnica ~ log(Stems) * Species + log(Stems) * Group + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m8)
check_overdispersion(m8)
#no overdispersion, interaction log(Stems) and group significant

m9 <- glmmTMB(Nr_Arnica ~ Stems * Species + Group + offset(log(nPoll))
              + (1|Site), family = nbinom2,
              data = comb_all2,
              control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m9)
check_overdispersion(m9)
#no overdispersion

#try with nPoll as predictor, not as offset
m10 <- glmmTMB(Nr_Arnica ~ (log(Stems) * Species + log(Stems) * Group) * log(nPoll)
               + (1|Site), family = nbinom2,
               data = comb_all2,
               control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m10)
check_overdispersion(m10)


mlist = list(m4, m5, m6, m7, m8, m9, m10)
AICTab = AIC(m4, m5, m6, m7, m8, m9, m10) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#ranked m10 < m8 < m9 < m7 < m5 < m4 < m6.  
#Negative binomial over poisson over binomial.

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
r.squaredGLMM(m5)
r.squaredGLMM(m6)
r.squaredGLMM(m7)
r.squaredGLMM(m8)
r.squaredGLMM(m9)
r.squaredGLMM(m10)
#does not currently work for m1, m4 and m6 slightly higher r2 than m5, then m7.
#negative binomial models much lower r2.

#model binomial----
m_species <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ log(Stems) * Species + log(Stems) * Group
                     + (1|Site), family = binomial,
                     data = comb_all2,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species)

eff_species <- effect("log(Stems)",m_species, xlevels = 50)  
eff.plot(eff_species, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "binomial model",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_species)
#no ovdispersion
#qqnorm(resid(m_species))
hist(resid(m_species))
#qqplot does not look too good, residual normality seems fine


residuals_species <- simulateResiduals(fittedModel = m_species)
plot(residuals_species)
testOutliers(residuals_species)
#all significant, model fit issues?

#model negative binomial----
m_species2 <- glmmTMB(Nr_Arnica ~ (log(Stems) * Species + log(Stems) * Group) * log(nPoll)
                      + (1|Site), family = nbinom2,
                      data = comb_all2,
                      control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species2)

eff_species2 <- effect(c("log(Stems)", "Group"),m_species2, xlevels = 50)
eff_species2 <- effect("log(Stems):Group",m_species2, xlevels = 50)
eff.plot(eff_species2, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "negative binomial model",
         ylim.data = T, overlay = F, 
         col.data = c(1,2))

#test if model assumptions are met and test model for fit:
check_overdispersion(m_species2)
#no ovdispersion
#qqnorm(resid(m_species2))
hist(resid(m_species2))
#qqplot and histogram of residuals look terrible

residuals_species2 <- simulateResiduals(fittedModel = m_species2)
plot(residuals_species2)
testOutliers(residuals_species2)
#outlier test non-significant, residuals vs predicted looks a liitle better than 
#for binomial model. KS test non-significant. Model fit issues?

#model poisson----
m_species3 <- glmmTMB(Nr_Arnica ~ log(Stems) * Species + log(Stems) * Group + offset(log(nPoll))
                      + (1|Site), family = poisson,
                      data = comb_all2,
                      control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))
summary(m_species3)

eff_species3 <- effect("log(Stems)",m_species3, xlevels = 50)  
eff.plot(eff_species3, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "poisson model",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_species3)
#overdispersion
#qqnorm(resid(m_species3))
hist(resid(m_species3))
#qqplot and histogram do not look good but better than negative binomial

residuals_species3 <- simulateResiduals(fittedModel = m_species3)
plot(residuals_species3)
testOutliers(residuals_species3)
#outlier test significant, residuals vs predicted looks worse than for negative
#binomial model. KS test significant. Model fit issues?

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


