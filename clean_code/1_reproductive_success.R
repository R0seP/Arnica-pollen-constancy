##This script explores the impact of the population size of Arnica on the reproductive
#success of the plant, proxied by the seed filling rate.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
seed_data <- read.csv("seed_data.csv", h = T)

#model selection----
#binomial models:
m1 <- glmmTMB(cbind(Filled, Not_Filled) ~ Stems + (1|Site),
              data = seed_data, family = binomial)

m2 <- glmmTMB(cbind(Filled, Not_Filled) ~ log(Stems) + (1|Site),
                      data = seed_data, family = binomial)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#m1 without log ranked higher than m2

check_overdispersion(m1)
check_overdispersion(m2)
#both overdispersed

#negative binomial models:
m3 <- glmmTMB(Filled ~ Stems +offset(log(Total_seeds))  + (1|Site),
              data = seed_data, family = nbinom2)

m4 <- glmmTMB(Filled ~ log(Stems) +offset(log(Total_seeds))  + (1|Site),
              data = seed_data, family = nbinom2)

m5 <- glmmTMB(Filled ~ log(Stems) + log(Total_seeds)  + (1|Site),
              data = seed_data, family = nbinom2)


mlist = list(m3, m4, m5)
AICTab = AIC(m3, m4, m5) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m4 with log ranked highest

check_overdispersion(m3)
check_overdispersion(m4)
check_overdispersion(m5)
#no overdispersion

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
r.squaredGLMM(m5)
#model m5 with nPoll as predictor and no offset by far highest r2

#binomial model----
#modelling the reproductive success of Arnica montana with a binomial model, 
#filled seeds are counted as "successes", empty seeds as "failures"
m_success <- glmmTMB(cbind(Filled, Not_Filled) ~ Stems + (1|Site),
                     data = seed_data, family = binomial)

summary(m_success)

#r-squared value:
r.squaredGLMM(m_success)

#test model fit:
eff_success <- effect("Stems",m_success, xlevels = 50)  
eff.plot(eff_success, plotdata = T,
         ylab = "Proportion of filled seeds",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


#test if model assumptions are met and test model for fit:
check_overdispersion(m_success)
#model overdispersed

hist(resid(m_success))

residuals_success <- simulateResiduals(fittedModel = m_success)
plot(residuals_success)
testOutliers(residuals_success)
#KS test significant, significant res vs pred, and significant outliers (p = 0.08)
#model fit issues!


#n.b. model proportion----
#fit a negative binomial model which uses the proportion of filled seeds as a response

##dummy model----
#For some reason, effect() cannot handle negative binomial models with only one predictor
#(other model type with only one predictor work fine, as well as nb models with several predictors).

#Therefore, I introduce a dummy parameter "Nonsense" to the model.

#"Nonsense" is uniform and does not affect the other model parameters, but it fixes the bug and 
#lets me plot my effect plot.

seed_data_small <- seed_data[,c(1,2,5,8)] #create smaller data frame with only necessary columns
seed_data_small$Nonsense <- rep(1,381) #add a second, uniform predictor

m_success2 <- glmmTMB(Filled ~ log(Stems) + Nonsense + offset(log(Total_seeds))  + (1|Site),
                      data = seed_data_small, family = nbinom2)

eff_success2 <- effect("log(Stems)", m_success2, xlevels = 50)
eff.plot(eff_success2, plotdata = T,
         ylab = "Number of filled seeds",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)



##model (this model is used in thesis!)----

m_success2 <- glmmTMB(Filled ~ log(Stems) + offset(log(Total_seeds))  + (1|Site),
              data = seed_data, family = nbinom2)
summary(m_success2)

#r-squared value:
r.squaredGLMM(m_success2)


#test if model assumptions are met and test model for fit:
check_overdispersion(m_success2)
#model no longer overdispersed!

hist(resid(m_success2))

residuals_success2 <- simulateResiduals(fittedModel = m_success2)
plot(residuals_success2)
testOutliers(residuals_success2)
#significant KS deviation and significant quantile deviations, but no longer any
#significant outliers -> better than binomial model?


##alternative model----
#use nPoll as covariate not as offset
m_success3 <- glmmTMB(Filled ~ log(Stems) + log(Total_seeds)  + (1|Site),
                      data = seed_data, family = nbinom2)

summary(m_success3)
#stems still significantly positive, nPoll also significant

#test if model assumptions are met and test model for fit:
check_overdispersion(m_success3)

hist(resid(m_success3))

residuals_success3 <- simulateResiduals(fittedModel = m_success3)
plot(residuals_success3)
testOutliers(residuals_success3)
#KS and res. vs pred. significant, looks very similar to m_success2 -> no real 
#improvement, return to previous model


##effect sizes----
#predict proportion of filled seeds for 10, 100, and 500 stems
pred_data_success2 <- data.frame(Stems = c(10, 100, 500), Total_seeds = mean(seed_data$Total_seeds))
predictions_success2 <- predict(m_success2, newdata = pred_data_success2, type = "response", re.form = NA, se.fit = T)

pred_data_success2$Pred.Filled <- predictions_success2$fit
pred_data_success2$Pred.Filled_SE <- predictions_success2$se.fit
pred_data_success2
#for 10 stems, 59.4 (+/- 6.1) seeds are expected to be filled, for 100 stems
#77.2 (+/- 4.4) seeds, and for 500 stems 92.7 (+/- 9.5) 
#seeds are expected to be filled per flower head, at an average of 105 seeds per
#flower head and averaged over all sites.

#calculate proportions:
pred_data_success2_prop <- data.frame(Stems = c(10, 100, 500), Total_seeds = mean(seed_data$Total_seeds))
pred_data_success2_prop$Mean_prop <- pred_data_success2$Pred.Filled/mean(seed_data$Total_seeds)
pred_data_success2_prop$Lower_mean <- (pred_data_success2$Pred.Filled - pred_data_success2$Pred.Filled_SE)/mean(seed_data$Total_seeds)
pred_data_success2_prop$Upper_mean <- (pred_data_success2$Pred.Filled + pred_data_success2$Pred.Filled_SE)/mean(seed_data$Total_seeds)
pred_data_success2_prop$Mean_prop_SE <- pred_data_success2_prop$Upper_mean - pred_data_success2_prop$Lower_mean
pred_data_success2_prop
#for 10 stems, 56.7 (+/- 11.6)% of seeds are expected to be filled, for 100 stems
#73.7 (+/- 8.4)% of seeds, and for 500 stems 88.6 (+/- 18.1)% of 
#seeds are expected to be filled per flower head, at an average of 105 seeds per
#flower head and averaged over all sites.

##additional plot----
seed_data_small$Stems <- as.factor(seed_data_small$Stems)
seed_data_small$Filling_Ratio <- seed_data$Filling_Ratio

#data only of proportion of filled seeds (like binomial model plot but without model)
par(bty = "l")
plot(seed_data$Filling_Ratio ~ seed_data$Stems,
     xlab = "Population size Arnica (Nr Stems)",
     ylab = "Proportion of filled seeds",
     col = rgb(.15,.15,.95,.25),
     pch = 21, bg = rgb(.15,.15,.95,.25))



#n.b. model number----
#just consider number of filled seeds, regardless of total number of seeds

seed_data_small$Stems <- as.numeric(seed_data_small$Stems)

##dummy model----
m_success4 <- glmmTMB(Filled ~ log(Stems)  + Nonsense + (1|Site),
                      data = seed_data_small, family = nbinom1)
#add nonsense parameter again for later plotting

summary(m_success4)
#stems still significantly positive, nPoll also significant

eff_success4 <- effect("log(Stems)", m_success4, xlevels = 50)
eff.plot(eff_success4, plotdata = T,
         ylab = "Number of filled seeds",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


##model (this model is used in thesis!)----
m_success4 <- glmmTMB(Filled ~ log(Stems)  + (1|Site),
                      data = seed_data, family = nbinom1)

summary(m_success4)

#r-squared value
r.squaredGLMM(m_success4)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_success4)
#was underdispersed, therefore switched to nbinom1

hist(resid(m_success4))

residuals_success4 <- simulateResiduals(fittedModel = m_success4)
plot(residuals_success4)
testOutliers(residuals_success4)
#KS and res. vs pred. significant, might look better than for other models?

##effect sizes----
pred_data_success4 <- data.frame(Stems = c(10, 100, 500))
predictions_success4 <- predict(m_success4, newdata = pred_data_success4, type = "response", re.form = NA, se.fit = T)

pred_data_success4$Pred.Filled <- predictions_success4$fit
pred_data_success4$Pred.Filled_SE <- predictions_success4$se.fit
pred_data_success4
#for 10 stems, 52.5 (+/- 6.5) seeds are expected to be filled, for 100 stems
#77.3 (+/- 5.3) seeds, and for 500 stems 101.4 (+/- 12.4) 
#seeds are expected to be filled per flower head.

