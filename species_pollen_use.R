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

#data visualization----
#calculate the mean P_ASTE.Arnica.montana for each species
species_means <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE)) %>%
  arrange(desc(mean_value))

#reorder the Species factor levels based on the calculated means
comb_all2$Species <- factor(comb_all2$Species, levels = species_means$Species)

#boxplot of percentage of Arnica carried
bp1 <- ggplot(comb_all2, aes(fill=Group, y=P_ASTE.Arnica.montana, x=Species)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(legend.text = element_text(size = 12), # increase legend text size
        legend.title = element_text(size = 14), # increase legend title size
        axis.text = element_text(size = 12), # increase axis text size
        axis.title = element_text(size = 14)) + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Species", y = "% Arnica montana pollen carried")

#calculate the mean Nr_Arnica for each species
species_means <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(Nr_Arnica, na.rm = TRUE)) %>%
  arrange(desc(mean_value))

#reorder the Species factor levels based on the calculated means
comb_all2$Species <- factor(comb_all2$Species, levels = species_means$Species)

#boxplot of number of Arnica pollen carried
bp2 <- ggplot(comb_all2, aes(fill=Group, y=Nr_Arnica, x=Species)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(legend.text = element_text(size = 12), # increase legend text size
        legend.title = element_text(size = 14), # increase legend title size
        axis.text = element_text(size = 12), # increase axis text size
        axis.title = element_text(size = 14)) + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Species", y = "Nr of Arnica montana pollen carried")

library(gridExtra)
grid.arrange(bp1, bp2, ncol = 1)


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


#model negative binomial----
#(m9.3 from model selection)
#relevel species factor so that Stenurella melanura is used as baseline (average intercept)
comb_all2$Species <- relevel(comb_all2$Species, ref = "Stenurella melanura")
#set the species with mean number of  log(Arnica carried) as reference 

m_species <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                     data = comb_all2, family = nbinom1,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))

summary(m_species)

check_overdispersion(m_species) #no overdispersion
check_collinearity(m_species) #no highly correlated predictors

eff_species <- effect(c("Stems"),m_species, xlevels = 50)
eff.plot(eff_species, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, 
         col.data = 3)

#find number of degrees of freedom
library(car)
anova <- Anova(m_species, type = "III")  # Type III ANOVA table
print(anova)

#test if model assumptions are met and test model for fit:


#qqnorm(resid(m_species2))
hist(resid(m_species)) #distribution of residuals high peak and more negative

residuals_species <- simulateResiduals(fittedModel = m_species)
plot(residuals_species)
testOutliers(residuals_species)


#effect sizes----
#r-squared value:
r.squaredGLMM(m_species)

#predict the percentage of Arnica pollen for species in different groups for 10, 100, 500 stems
pred_data_species <- data.frame(Stems = c(rep(10,44),rep(100,44),rep(500,44)), 
                                Species = rep(imp_species$Species,6), 
                                Group = c(rep("area",22), rep("flower",22)),
                                nPoll = rep(mean(comb_all2$nPoll), 132))
predictions_species <- predict(m_species, newdata = pred_data_species, type = "response", re.form = NA, se.fit = T)

pred_data_species$Pred.Arnica <- predictions_species$fit
pred_data_species$Pred.Arnica_SE <- predictions_species$se.fit
pred_data_species


#models by species----
###Andrena sp----
Andrena <- subset(comb_all2, Species == "Andrena sp")

mAndrena <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                    + (1|Site), family = nbinom1,
                    data = Andrena)

check_overdispersion(mAndrena)
summary(mAndrena)
r.squaredGLMM(mAndrena)

eff_Andrena <- effect(c("Stems"),mAndrena, xlevels = 50)
eff.plot(eff_Andrena, plotdata = T,
         ylab = "Nr of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Adrena sp",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Andrena <- simulateResiduals(fittedModel = mAndrena)
plot(residuals_Andrena)
testOutliers(residuals_Andrena)
#no fit issues


###Apis mellifera----
Apis <- subset(comb_all2, Species == "Apis mellifera")

mApis <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                  + (1|Site), family = nbinom1,
                  data = Apis)

check_overdispersion(mApis)
summary(mApis)
r.squaredGLMM(mApis)

eff_Apis <- effect(c("Stems"),mApis, xlevels = 50)
eff.plot(eff_Apis, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Apis mellifera",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Apis <- simulateResiduals(fittedModel = mApis)
plot(residuals_Apis)
testOutliers(residuals_Apis)
#significant KS and significant res vs pred, binomial model much better fit here

###Bombus pascuorum----
Bombus_pascuorum <- subset(comb_all2, Species == "Bombus pascuorum")

mBombus_pascuorum <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                             + (1|Site), family = nbinom1,
                             data = Bombus_pascuorum)

check_overdispersion(mBombus_pascuorum)
summary(mBombus_pascuorum)
r.squaredGLMM(mBombus_pascuorum)

eff_Bombus_pascuorum <- effect("Stems",mBombus_pascuorum, xlevels = 50)
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

mBombus_ruderarius <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                              + (1|Site), family = nbinom1,
                              data = Bombus_ruderarius)

check_overdispersion(mBombus_ruderarius)
summary(mBombus_ruderarius)
r.squaredGLMM(mBombus_ruderarius)

eff_Bombus_ruderarius <- effect("Stems",mBombus_ruderarius, xlevels = 50)
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

mBombus_terrestris <-glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                             + (1|Site), family = nbinom1,
                             data = Bombus_terrestris)

check_overdispersion(mBombus_terrestris)
summary(mBombus_terrestris)
r.squaredGLMM(mBombus_terrestris)

eff_Bombus_terrestris <- effect("Stems",mBombus_terrestris, xlevels = 50)
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

mCoenonympha <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                        + (1|Site), family = nbinom1,
                        data = Coenonympha)

check_overdispersion(mCoenonympha)
summary(mCoenonympha)
r.squaredGLMM(mCoenonympha)

eff_Coenonympha <- effect("Stems",mCoenonympha, xlevels = 50)
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

mDasytes <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                    + (1|Site), family = nbinom1,
                    data = Dasytes)

check_overdispersion(mDasytes)
summary(mDasytes)
r.squaredGLMM(mDasytes)

eff_Dasytes <- effect("Stems",mDasytes, xlevels = 50)
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

mEmpis_livida <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                         + (1|Site), family = nbinom1,
                         data = Empis_livida)

check_overdispersion(mEmpis_livida)
summary(mEmpis_livida)
r.squaredGLMM(mEmpis_livida)

eff_Empis_livida <- effect("Stems",mEmpis_livida, xlevels = 50)
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

mEmpis_tessellata <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                             + (1|Site), family = nbinom1,
                             data = Empis_tessellata)

check_overdispersion(mEmpis_tessellata)
summary(mEmpis_tessellata)
r.squaredGLMM(mEmpis_tessellata)

eff_Empis_tessellata <- effect("Stems",mEmpis_tessellata, xlevels = 50)
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

mEristalis <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                              + (1|Site), family = nbinom1,
                             data = Eristalis)

check_overdispersion(mEristalis)
summary(mEristalis)
r.squaredGLMM(mEristalis)

eff_Eristalis <- effect("Stems",mEristalis, xlevels = 50)
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

#negative binomial model
mEupeodes1 <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                     + (1|Site), family = nbinom1,
                     data = Eupeodes)

check_overdispersion(mEupeodes1)
summary(mEupeodes1)
r.squaredGLMM(mEupeodes1)

eff_Eupeodes1 <- effect("Stems",mEupeodes1, xlevels = 50)
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


###Helophilus pendulus----
Helophilus <- subset(comb_all2, Species == "Helophilus pendulus")

mHelophilus <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                       + (1|Site), family = nbinom1,
                       data = Helophilus)

check_overdispersion(mHelophilus)
summary(mHelophilus)
r.squaredGLMM(mHelophilus)

eff_Helophilus <- effect("Stems",mHelophilus, xlevels = 50)
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

mLasioglossum <-glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                        + (1|Site), family = nbinom1,
                        data = Lasioglossum)

check_overdispersion(mLasioglossum)
summary(mLasioglossum)
r.squaredGLMM(mLasioglossum)

eff_Lasioglossum <- effect("Stems",mLasioglossum, xlevels = 50)
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

mManiola <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                    + (1|Site), family = nbinom1,
                    data = Maniola)

check_overdispersion(mManiola)
summary(mManiola)
r.squaredGLMM(mManiola)

eff_Maniola <- effect("Stems",mManiola, xlevels = 50)
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

mMeligethes <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                       + (1|Site), family = nbinom1,
                       data = Meligethes)

check_overdispersion(mMeligethes)
summary(mMeligethes)
r.squaredGLMM(mMeligethes)

eff_Meligethes <- effect("Stems",mMeligethes, xlevels = 50)
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

mMerodon <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                    + (1|Site), family = nbinom1,
                    data = Merodon)

check_overdispersion(mMerodon)
summary(mMerodon)
r.squaredGLMM(mMerodon)

eff_Merodon <- effect("Stems",mMerodon, xlevels = 50)
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

mNomada <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                   + (1|Site), family = nbinom1,
                   data = Nomada)

check_overdispersion(mNomada)
summary(mNomada)
r.squaredGLMM(mNomada)

eff_Nomada <- effect("Stems",mNomada, xlevels = 50)
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

mOchlodes <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll))
                     + (1|Site), family = nbinom1,
                     data = Ochlodes)

check_overdispersion(mOchlodes)
summary(mOchlodes)
r.squaredGLMM(mOchlodes)

eff_Ochlodes <- effect("Stems",mOchlodes, xlevels = 50)
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

mOedemera <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                     + (1|Site), family = nbinom1,
                     data = Oedemera)

check_overdispersion(mOedemera)
summary(mOedemera)
r.squaredGLMM(mOedemera)

eff_Oedemera <- effect("Stems",mOedemera, xlevels = 50)
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

mPhyllopertha <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                         + (1|Site), family = nbinom1,
                         data = Phyllopertha)

check_overdispersion(mPhyllopertha)
summary(mPhyllopertha)
r.squaredGLMM(mPhyllopertha)

eff_Phyllopertha <- effect("Stems",mPhyllopertha, xlevels = 50)
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

mSphaerophoria <- glmmTMB(Nr_Arnica ~ Stems + Group + offset(log(nPoll))
                          + (1|Site), family = nbinom1,
                          data = Sphaerophoria)

check_overdispersion(mSphaerophoria)
summary(mSphaerophoria)
r.squaredGLMM(mSphaerophoria)

eff_Sphaerophoria <- effect("Stems",mSphaerophoria, xlevels = 50)
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

mStenurella <- glmmTMB(Nr_Arnica ~ Stems * Group + offset(log(nPoll))
                       + (1|Site), family = nbinom1,
                       data = Stenurella)

check_overdispersion(mStenurella)
summary(mStenurella)
r.squaredGLMM(mStenurella)

eff_Stenurella <- effect("Stems",mStenurella, xlevels = 50)
eff.plot(eff_Stenurella, plotdata = T,
         ylab = "Proportion of Arnica pollen carried",
         xlab = "Population size Arnica (Nr Stems)",
         main = "Stenurella mealnura",
         ylim.data = T, overlay = F, 
         col.data = 3)

residuals_Stenurella <- simulateResiduals(fittedModel = mStenurella)
plot(residuals_Stenurella)
testOutliers(residuals_Stenurella)


#model visualization----
comb_all2 <- na.omit(comb_all2)
comb_all2$predictions <- predict(m_species, type = "response")

theta_value <- sigma(m_species)
library(MASS)

ggplot(comb_all2, aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "", x = "  Nr Stems", y = "Predicted Nr Arnica pollen carried", color = "Group")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p0 <- ggplot(comb_all2, aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "All species", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

#by species:
p1 <- ggplot(comb_all2[comb_all2$Species == "Andrena sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Andrena sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p2 <- ggplot(comb_all2[comb_all2$Species == "Apis mellifera",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Apis mellifera", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p3 <- ggplot(comb_all2[comb_all2$Species == "Bombus pascuorum",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Bombus pascuorum", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p4 <- ggplot(comb_all2[comb_all2$Species == "Bombus ruderarius",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Bombus ruderarius", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p5 <- ggplot(comb_all2[comb_all2$Species == "Bombus terrestris",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Bombus terrestris", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p6 <- ggplot(comb_all2[comb_all2$Species == "Coenonympha pamphilus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Coenonympha pamphilus", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p7 <- ggplot(comb_all2[comb_all2$Species == "Dasytes niger",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Dasytes niger", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p8 <- ggplot(comb_all2[comb_all2$Species == "Empis livida",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Empis livida", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p9 <- ggplot(comb_all2[comb_all2$Species == "Empis tessellata",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Empis tesellata", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p10 <- ggplot(comb_all2[comb_all2$Species == "Eristalis sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Eristalis sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p11 <- ggplot(comb_all2[comb_all2$Species == "Eupeodes corollae",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Eupeodes corollae", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p12 <- ggplot(comb_all2[comb_all2$Species == "Helophilus pendulus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Helophilus pendulus", x = "", y = "Predicted Nr of Arnica pollen", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p13 <- ggplot(comb_all2[comb_all2$Species == "Lasioglossum sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Lasioglossum sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p14 <- ggplot(comb_all2[comb_all2$Species == "Maniola jurtina",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Maniola jurtina", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p15 <- ggplot(comb_all2[comb_all2$Species == "Meligethes sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Meligethes sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p16 <- ggplot(comb_all2[comb_all2$Species == "Merodon equestris",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Merodon equestris", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p17 <- ggplot(comb_all2[comb_all2$Species == "Nomada sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Nomada sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p18 <- ggplot(comb_all2[comb_all2$Species == "Ochlodes sylvanus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Ochlodes sylvanus", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p19 <- ggplot(comb_all2[comb_all2$Species == "Oedemera sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Oedemera sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p20 <- ggplot(comb_all2[comb_all2$Species == "Phyllopertha horticola",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Phyllopertha horticola", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p21 <- ggplot(comb_all2[comb_all2$Species == "Sphaerophoria sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Sphaerophoria sp", x = "Nr Stems", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p22 <- ggplot(comb_all2[comb_all2$Species == "Stenurella melanura",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "glm", method.args = list(family = negative.binomial(theta = theta_value))) +
  labs(title = "Stenurella melanura", x = "", y = "", color = "Group") +
  theme(legend.position = c(1.5, 0.5), legend.justification = c(1, 0.5)) +
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

grid.arrange(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
             p15, p16, p17, p18, p19, p20, p21, p22, ncol = 6)

