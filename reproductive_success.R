##This script explores the impact of the population size of Arnica on the reproductive
#success of the plant, proxied by the seed filling rate.

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


mlist = list(m3, m4)
AICTab = AIC(m3, m4) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m4 with log ranked higher

check_overdispersion(m3)
check_overdispersion(m4)
#no overdispersion

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
#model m1, binomial without log highest r2, which model to use?

#binomial model----
#modelling the reproductive success of Arnica montana with a binomial model, 
#filled seeds are counted as "successes", empty seeds as "failures"
m_success <- glmmTMB(cbind(Filled, Not_Filled) ~ Stems + (1|Site),
                     data = seed_data, family = binomial)

summary(m_success)

eff_success <- effect("Stems",m_success, xlevels = 50)  
eff.plot(eff_success, plotdata = T,
         ylab = "Proportion of filled seeds",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


#test if model assumptions are met and test model for fit:
check_overdispersion(m_success)
#model overdispersed
#qqnorm(resid(m_success))
hist(resid(m_success))
#qqplot and normality of residuals look good to me

residuals_success <- simulateResiduals(fittedModel = m_success)
plot(residuals_success)
testOutliers(residuals_success)
#KS test significant, significant res vs pred, and significant outliers (p = 0.08)
#there seem to be issues with the model, what can be done about it?

#effect sizes binomial----

#r-squared value:
r.squaredGLMM(m_success)

#predict proportion of filled seeds for 10, 100, and 500 stems
pred_data_success <- data.frame(Stems = c(10, 100, 500))
predictions_success <- predict(m_success, newdata = pred_data_success, type = "response", re.form = NA, se.fit = T)

pred_data_success$Pred.Filled <- predictions_success$fit
pred_data_success$Pred.Filled_SE <- predictions_success$se.fit
pred_data_success
#for 10 stems, 67.2 (+/- 4.6)% of seeds are expected to be filled, for 100 stems
#75.7 (+/- 3.2)% of seeds, and for 500 stems 91.8 (+/- 3.1)% of 
#seeds are expected to be filled per flower head, averaged over all sites.


#negative binomial model----
#seed_data$log_Stems <- log(seed_data$Stems)

#figure out why effect doesn't work
columns_with_zero_or_below <- which(apply(seed_data, 2, function(x) any(x <= 0)))
print(columns_with_zero_or_below)

rows_with_zero_or_below <- which(apply(seed_data, 1, function(row) any(row <= 0)))
print(rows_with_zero_or_below)

#entries_zero_or_below <- which(seed_data <= 0, arr.ind = TRUE)
#print(entries_zero_or_below)

seed_data_filtered <- seed_data[apply(seed_data, 1, function(row) all(row > 0)), ]
#returns no data

seed_data_filtered <- seed_data
seed_data_filtered$Stems <- pmax(pmin(seed_data_filtered$Stems, 
                                      max(seed_data_filtered$Stems)), 1)
seed_data_filtered$Head <- pmax(pmin(seed_data_filtered$Head, 
                                     max(seed_data_filtered$Head)), 1)
seed_data_filtered$Filled <- pmax(pmin(seed_data_filtered$Filled, 
                                     max(seed_data_filtered$Filled)), 1)
seed_data_filtered$Not_Filled <- pmax(pmin(seed_data_filtered$Not_Filled, 
                                     max(seed_data_filtered$Not_Filled)), 1)
seed_data_filtered$Total_seeds <- pmax(pmin(seed_data_filtered$Total_seeds, 
                                           max(seed_data_filtered$Total_seeds)), 1)

seed_data$Stems <- as.numeric(seed_data$Stems)
seed_data$Head <- as.numeric(seed_data$Head)
seed_data$Filled <- as.numeric(seed_data$Filled)
seed_data$Not_Filled <- as.numeric(seed_data$Not_Filled)
seed_data$Total_seeds <- as.numeric(seed_data$Total_seeds)
seed_data$Site <- as.factor(seed_data$Site)
str(seed_data)

m_success2 <- glmmTMB(Filled ~ log(Stems) + offset(log(Total_seeds))  + (1|Site),
              data = seed_data, family = nbinom2)
summary(m_success2)

eff_success2 <- effect("log(Stems)", m_success2, xlevels = 50)
eff.plot(effect_plot, plotdata = T,
         ylab = "Number of filled seeds",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
check_overdispersion(m_success2)
#model no longer overdispersed!
#qqnorm(resid(m_success2))
hist(resid(m_success2))
#qqplot and normality of residuals look good to me, maybe residuals skewed 
#towards lower values

residuals_success2 <- simulateResiduals(fittedModel = m_success2)
plot(residuals_success2)
testOutliers(residuals_success2)
#significant KS deviation and significant quantile deviations, but no longer any
#significant outliers -> better than binomial model?
#Just because of big sample size and model family?

#effect sizes n.b.----

#r-squared value:
r.squaredGLMM(m_success2)

#predict number of filled seeds for 10, 100, and 500 stems
pred_data_success2 <- data.frame(Stems = c(10, 100, 500), Total_seeds = mean(seed_data$Total_seeds))
predictions_success2 <- predict(m_success2, newdata = pred_data_success2, type = "response", re.form = NA, se.fit = T)

pred_data_success2$Pred.Filled <- predictions_success2$fit
pred_data_success2$Pred.Filled_SE <- predictions_success2$se.fit
pred_data_success2
#for 10 stems, 59.4 (+/- 6.1)% of seeds are expected to be filled, for 100 stems
#77.2 (+/- 4.4)% of seeds, and for 500 stems 92.7 (+/- 9.5)% of 
#seeds are expected to be filled per flower head, at an average of 105 seeds per
#flower head and averaged over all sites.


#alternative way of creating a plot since effect() currently doesn't work, with fixed offset term
library(ggeffects)
effect_plot <- ggpredict(m_success2, terms = "Stems", condition = c(Total_seeds = 1))
plot(effect_plot)
