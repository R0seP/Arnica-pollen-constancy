#This script investigates whether the Bray-Curtis dissimilarity between pollinators
#caught on Arnica flowers and pollinators caught in the area at a site is affected
#by the population size of Arncia montana at that site. 

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)

library(vegan)
library(tidyverse)

#model all species----
###data preparation----
data <- read.csv("species_numbers_persite.csv")

#split data by site
sites <- split(data, data$Site)

#calculate Bray-Curtis dissimilarity for each site
bray_curtis <- sapply(sites, function(site_data) {
  vegdist(site_data[, -c(1, 2)], method = "bray")
})

#extract dissimilarity between the two groups at each site
bray_curtis_values <- sapply(bray_curtis, function(x) x[1])
bray_curtis_values

#add Bray-Curtis dissimilarity to other data
comb1$Bray_dissim <- bray_curtis_values

#clean up and sort comb1
comb1 <- comb1[,-3]
comb1 <- comb1[,c(1:2,30,3:29)]

###model selection----
m1 <- glmmTMB(Bray_dissim ~ Stems, data = comb1)
m2 <- glmmTMB(Bray_dissim ~ log(Stems), data = comb1)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#m1 without log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#m1 slightly higher r2, both very low r2

###model----
m_bray <- glmmTMB(Bray_dissim ~ Stems, data = comb1)

summary(m_bray)

eff_bray <- effect("Stems",m_bray, xlevels = 50)  
eff.plot(eff_bray, plotdata = T,
         ylab = "Bray-Curtis dissimilarity",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_bray)) #very jagged
hist(resid(m_bray)) #does not look normally distributed

residuals_bray <- simulateResiduals(fittedModel = m_bray)
plot(residuals_bray)
testOutliers(residuals_bray)
#nothing significant, fit okay and deviations in visual inspection because of
#few data points?


###effect sizes----
#r-squared value:
r.squaredGLMM(m_bray)


#predict Bray-Curtis distance between pollinators caught on Arnica flowers and pollinators caught in the area
#for 10, 100, and 500 stems
pred_data_bray <- data.frame(Stems = c(10, 100, 500))
predictions_bray <- predict(m_bray, newdata = pred_data_bray, type = "response", se.fit = TRUE)

pred_data_bray$Pred.Dissim <- predictions_bray$fit
pred_data_bray$Pred.Dissim_SE <- predictions_bray$se.fit
pred_data_bray
#for 10 stems, a dissimilarity of 85.7 (+/- 3.1)% is expected, for 100 stems
#a dissimilarity of 84.8 (+/- 2.4)%, and for 500 stems a dissimilarity of 
#80.9 (+/- 6.1)% is expected.


#model imp species only----
###data preparation----
data2 <- read.csv("species_numbers_persite.csv")
data2 <- data2[, c(1,2,8,10,16,18,21,31,35,38,39,42,45,46,50,53,59,63,64,66,
                 76,78)]
#remove all but species with 5 or more observations
#remove Nomada because of no association
#remove Meligethes because of potential sampling bias
#data set now comparable to "Arnica associated" and "Area associated" analysis

#split data by site
sites2 <- split(data2, data2$Site)

#calculate Bray-Curtis dissimilarity for each site
bray_curtis <- sapply(sites2, function(site_data) {
  vegdist(site_data[, -c(1, 2)], method = "bray")
})

#extract dissimilarity between the two groups at each site
bray_curtis_values <- sapply(bray_curtis, function(x) x[1])
bray_curtis_values

#add Bray-Curtis dissimilarity to other data
comb1$Bray_dissim2 <- bray_curtis_values

#clean up and sort comb1
comb1 <- comb1[,c(1:3,31,4:30)]

###model selection----
m1 <- glmmTMB(Bray_dissim2 ~ Stems, data = comb1)
m2 <- glmmTMB(Bray_dissim2 ~ log(Stems), data = comb1)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#m1 without log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#m1 slightly higher r2, both very low r2

###model----
m_bray2 <- glmmTMB(Bray_dissim2 ~ Stems, data = comb1)

summary(m_bray2)

eff_bray2 <- effect("Stems",m_bray2, xlevels = 50)  
eff.plot(eff_bray2, plotdata = T,
         ylab = "Bray-Curtis dissimilarity",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_bray2)) #loks fine
hist(resid(m_bray2)) #does not look normally distributed

residuals_bray2 <- simulateResiduals(fittedModel = m_bray2)
plot(residuals_bray2)
testOutliers(residuals_bray2)
#nothing significant, fit okay and deviations in visual inspection because of
#few data points?


###effect sizes----
#r-squared value:
r.squaredGLMM(m_bray2)


#predict Bray-Curtis dissimilarity only with imp species
#for 10, 100, and 500 stems
pred_data_bray2 <- data.frame(Stems = c(10, 100, 500))
predictions_bray2 <- predict(m_bray2, newdata = pred_data_bray2, type = "response", se.fit = TRUE)

pred_data_bray2$Pred.Dissim <- predictions_bray2$fit
pred_data_bray2$Pred.Dissim_SE <- predictions_bray2$se.fit
pred_data_bray2
#for 10 stems, a dissimilarity of 83.1 (+/- 3.5)% is expected, for 100 stems
#a dissimilarity of 81.7 (+/- 2.7)%, and for 500 stems a dissimilarity of 
#75.3 (+/- 7.0)% is expected.

