#This script looks at the effect of population size of Arnica on the number of 
#pollinators caught while controlling for environmental and temporal effects.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
poll_nr <- read.csv("poll_nr.csv", h = T)

poll_nr$Date_from_start_factor <- as.factor(poll_nr$Date_from_start)
poll_nr$Poll_per_hr <- round(poll_nr$Poll_per_hr)


#model----
m_poll_nr <- glmmTMB(Poll_per_hr ~ Stems + Date_from_start
                     + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)

summary(m_poll_nr)
 
#r-squared:
r.squaredGLMM(m_poll_nr)


eff_poll_nr <- effect("Stems",m_poll_nr, xlevels = 50)  
eff.plot(eff_poll_nr, plotdata = T,
         ylab = "Number of pollinators per hour",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

eff_poll_nr2 <- effect("Date_from_start",m_poll_nr, xlevels = 50)  
eff.plot(eff_poll_nr2, plotdata = T,
         ylab = "Number of pollinators per hour",
         xlab = "Sampling day from start",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#model test----
#test if model assumptions are met and test model for fit:
check_overdispersion(m_poll_nr)

hist(resid(m_poll_nr))

ks.test(resid(m_poll_nr), "pnorm", mean = mean(resid(m_poll_nr)), sd = sd(resid(m_poll_nr)))
#ks test non-significant

residuals_poll_nr <- simulateResiduals(fittedModel = m_poll_nr)
plot(residuals_poll_nr)
testOutliers(residuals_poll_nr)
#no outliers significant, but combined adjusted quantile test significant, 
#there are residuals vs predictions quantile deviations


#effect size----

#predict numbers of pollinators caught per hour
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


#for sampling day 1, 10, and 20
pred_data_poll_nr2 <- data.frame(Stems = rep(mean(poll_nr$Stems),3), 
                                Date_from_start = c(1,10,20) )
predictions_poll_nr2 <- predict(m_poll_nr, newdata = pred_data_poll_nr2, type = "response", re.form = NA, se.fit = T)

pred_data_poll_nr2$Pred.Polli_Nr <- predictions_poll_nr2$fit
pred_data_poll_nr2$Pred.POlli_Nr_SE <- predictions_poll_nr2$se.fit
pred_data_poll_nr2
#at sampling day 1, the model predicts 8.9 (+/- 1.2) pollinators to be caught per hour, 
#at sampling day 10, it predicts 10.5 (+/- 0.9) pollinators/hr,
#at sampling day 20 12.5 (+/- 0.9) pollinators/hr. 