#This script extracts the number of pollinators caught in the area that also use 
#Arnica (defined by percentage of Arnica over 5%) and models whether that number is influenced
#by the number of Stems in the area.

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

#data preparation----
#create data frame with counts of area pollinators that carry >5% of Arnica
area_subset <- subset(comb_all, Group == "area")    #only area samples
Arnica5 <- subset(area_subset, P_ASTE.Arnica.montana >= 0.05)   #only samples with >= 5% Arnica
Arnica5_counts <- Arnica5 %>%
  count(Site)     #create df that contains counts of occurrences of area pollinators with >= 5% Arnica
missing_sites <- data.frame(Site = c(3,5,13,15), n = rep(0,4)) 
Arnica5_counts <- rbind(Arnica5_counts, missing_sites) #add sites that have no occurrences
Arnica5_counts <- Arnica5_counts[order(Arnica5_counts$Site), ] #sort by site
Arnica5_counts$Stems <- comb1$Stems #add Arnica population size

#model selection----
m1 <- glmmTMB(n ~ Stems, data = Arnica5_counts, family = poisson)
m2 <- glmmTMB(n ~ log(Stems), data = Arnica5_counts, family = poisson)

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
#model m1 without log also slightly higher r2, use that model!

#model----
m_Arnica5_area <- glmmTMB(n ~ Stems, data = Arnica5_counts, family = poisson)
summary(m_Arnica5_area)
check_overdispersion(m_Arnica5_area)
#no overdispersion

eff_Arnica5_area <- effect("Stems",m_Arnica5_area, xlevels = 50)  
eff.plot(eff_Arnica5_area, plotdata = T,
         ylab = "Number of pollinators caught in Area visiting Arnica",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#find number of degrees of freedom
library(car)
anova <- Anova(m_Arnica5_area, type = "III")  # Type III ANOVA table
print(anova)


#model test----
#test if model assumptions are met and test model for fit:
qqnorm(resid(m_Arnica5_area)) #issues at lower theoretical quantiles
hist(resid(m_Arnica5_area)) #residuals do not look normally distributed!

ks.test(resid(m_Arnica5_area), "pnorm", mean = mean(resid(m_Arnica5_area)), sd = sd(resid(m_Arnica5_area)))
#ks test non-significant, non-normality then fine?

residuals_Arnica5_area <- simulateResiduals(fittedModel = m_Arnica5_area)
plot(residuals_Arnica5_area)
testOutliers(residuals_Arnica5_area)
#no significant tests, model fine?

#effect size----
r.squaredGLMM(m_Arnica5_area)
#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_Arnica5_area <- data.frame(Stems = c(10, 100, 500))
predictions_Arnica5_area <- predict(m_Arnica5_area, newdata = pred_data_Arnica5_area, type = "response", re.form = NA, se.fit = T)

pred_data_Arnica5_area$Pred.Polli_Nr <- predictions_Arnica5_area$fit
pred_data_Arnica5_area$Pred.POlli_Nr_SE <- predictions_Arnica5_area$se.fit
pred_data_Arnica5_area
#for 10 stems, 2.1 (+/- 0.4) area-caught pollinators are expected to also have
#visited Arnica, for 100 stems 2.3 (+/- 0.4) pollinators, and for 500 stems
#4.2 (+/- 1.2) pollinators caught in the area are expected to have visited Arnica as well.


#alternative models----
Arnica5$Nr_Arnica <- round(Arnica5$Nr_Arnica)
m1 <- glmmTMB(Nr_Arnica ~ Stems + Species + offset(log(nPoll)), data = Arnica5,
              family = nbinom2)
summary(m1)
r.squaredGLMM(m1)

m2 <- glmmTMB(Nr_Arnica ~ Stems + offset(log(nPoll)), data = Arnica5,
              family = nbinom2)
summary(m2)
r.squaredGLMM(m2)
#model with species explains nearly all the variance in the data, while model 
#with only stems barely explains variance

eff_Arnica5_area2 <- effect("Stems",m1, xlevels = 50)  
eff.plot(eff_Arnica5_area2, plotdata = T,
         ylab = "Number of pollinators caught in Area visiting Arnica",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)
