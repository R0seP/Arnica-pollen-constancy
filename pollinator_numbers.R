#This script looks at the effect of population size of Arnica on the number of 
#pollinators caught while controlling for environmental and temporal effects.


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


#How is the number of pollinators caught per hour influenced by the population size of Arnica?

poll_nr$Date_from_start_factor <- as.factor(poll_nr$Date_from_start)
poll_nr$Poll_per_hr <- round(poll_nr$Poll_per_hr)

#model selection----
m0 <- glmmTMB(Poll_per_hr ~ 1, data = poll_nr, family = poisson)
m1 <- glmmTMB(Poll_per_hr ~ Stems * prec_mm.day * Date_from_start
                     + (1|time_quali) + (1|Site), data = poll_nr,
                     family = poisson)
m2 <- glmmTMB(Poll_per_hr ~ Stems + prec_mm.day * Date_from_start
              + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)
m3 <- glmmTMB(Poll_per_hr ~ Stems + prec_mm.day + Date_from_start
              + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)
m4 <- glmmTMB(Poll_per_hr ~ Stems * prec_mm.day + (1|Site)
              + (Date_from_start_factor|time_quali), data = poll_nr, family = poisson)
m5 <- glmmTMB(Poll_per_hr ~ Stems + prec_mm.day + (1|Site)
              + (Date_from_start_factor|time_quali), data = poll_nr, family = poisson)
m6 <- glmmTMB(Poll_per_hr ~ Stems + Date_from_start
              + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)
m7 <- glmmTMB(Poll_per_hr ~ Stems + prec_mm.day
              + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)
#models 1, 4 and 5 caused issues in the model selection and were ranked last,
#therefore they are removed from the selection process.
#also potential issues with treating Date_from_start as a factor!

mlist = list(m0, m2, m3, m6, m7)
AICTab = AIC(m0, m2, m3, m6, m7) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#m6 ranked highest, use m6!

#test whether m6 gets better if the log of Stems is taken:
m1 <- glmmTMB(Poll_per_hr ~ Stems + Date_from_start
             + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)
m2 <- glmmTMB(Poll_per_hr ~ log(Stems) + Date_from_start
              + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)

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
#model m1 without log also  a little bit higher than r2, use that model!


#model----
m_poll_nr <- glmmTMB(Poll_per_hr ~ Stems + Date_from_start
                     + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)

summary(m_poll_nr)
check_overdispersion(m_poll_nr)
#no overdispersion

eff_poll_nr <- effect("Stems",m_poll_nr, xlevels = 50)  
eff.plot(eff_poll_nr, plotdata = T,
         ylab = "Number of pollinators per hour",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


#model test----
#test if model assumptions are met and test model for fit:
qqnorm(resid(m_poll_nr)) #issues at higher theoretical quantiles
hist(resid(m_poll_nr)) #very normal residual distribution except for few high outliers

ks.test(resid(m_poll_nr), "pnorm", mean = mean(resid(m_poll_nr)), sd = sd(resid(m_poll_nr)))
#ks test non-significant, non-normality then fine?

residuals_poll_nr <- simulateResiduals(fittedModel = m_poll_nr)
plot(residuals_poll_nr)
testOutliers(residuals_poll_nr)
#no outliers significant, but combined adjusted quantile test significant, 
#there are residuals vs predictions quantile deviations


#effect size----
r.squaredGLMM(m_poll_nr)
#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_poll_nr <- data.frame(Stems = c(10, 100, 500), 
                                Date_from_start = rep(mean(poll_nr$Date_from_start),3))
predictions_poll_nr <- predict(m_poll_nr, newdata = pred_data_poll_nr, type = "response", re.form = NA, se.fit = T)

pred_data_poll_nr$Pred.Polli_Nr <- predictions_poll_nr$fit
pred_data_poll_nr$Pred.POlli_Nr_SE <- predictions_poll_nr$se.fit
pred_data_poll_nr
#for 10 stems, 11.4 (+/- 0.9) pollinators are expected to be caught per hour, 
#for 100 stems 11.7 (+/- 0.8) pollinators/hour, and for 500 stems
# 13.0 (+/- 1.8) pollinators per hour are expected to be caught,
#averaged over the number of days from start and the time of the day.
