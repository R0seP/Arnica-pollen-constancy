#This script calculates the effective number of pollen species for each sample 
#as Shannon diversity. It then models the effect of population size of Arnica 
#on that diversity.

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")

library(vegan)
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
adj <- read.csv("adj_new.csv", h = T)


#prepare data
pollen <- adj[,c(3,8:42)]
pollen <- na.omit(pollen)

#Shannon diversity----
#calculate Shannon (Renyi?) diversity of the pollen samples
distances <- renyi(pollen[,2:36], scales = c(0,1)) #calculate species richness and Shannon entropy (could do Shannon diversity directly with hill = T)
distances$filn <- pollen$filn
distances <- distances %>% select(,c("filn","0","1"))
distances <- distances %>%
  rename(Species_richness = "0", Shannon_entropy = "1")
distances$Shannon_diversity = 2^distances$Shannon_entropy #calculate Shannon diversity (as 2^H (H = Shannon entropy))
dist_data <- inner_join(distances, comb_all, by = "filn")


#model selection----
m1 <- glmmTMB(Shannon_diversity ~ Stems * Group + (1|Site) + (1|Species), 
              data = dist_data)
m2 <- glmmTMB(Shannon_diversity ~ log(Stems) * Group + (1|Site) + (1|Species), 
              data = dist_data)
m3 <- glmmTMB(Shannon_diversity ~ Stems + (1|Site) + (1|Species), 
              data = dist_data)
m4 <- glmmTMB(Shannon_diversity ~ log(Stems) + (1|Site) + (1|Species), 
              data = dist_data)
m5 <- glmmTMB(Shannon_diversity ~ log(Stems) + Group + (1|Site) + (1|Species), 
              data = dist_data)
m6 <- glmmTMB(Shannon_diversity ~ log(Stems) + I(log(Stems)^2) + (1|Site) + (1|Species), 
              data = dist_data)
m7 <- glmmTMB(Shannon_diversity ~ log(Stems) + I(log(Stems)^2) + Group 
              + (1|Site) + (1|Species), data = dist_data)

mlist = list(m1, m2, m3, m4, m5, m6, m7)
AICTab = AIC(m1, m2, m3, m4, m5, m6, m7) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#AIC m7 < m5 < m2 < m1 < m4 < m6 < m3

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
r.squaredGLMM(m5)
r.squaredGLMM(m6)
r.squaredGLMM(m7)
#r2: m7 > m2 > m5 > m1 > m6 > m4 > m3

#model----
#Model the effect of Arnica population size on the effective number of species a 
#pollinator was carrying

###non-quadratic----
#model ranked best without quaratic term (m5)
m_Shannon <- glmmTMB(Shannon_diversity ~ log(Stems) + Group + (1|Site) + (1|Species), 
                    data = dist_data)
summary(m_Shannon)

eff_Shannon <- effect("log(Stems)",m_Shannon, xlevels = 50)  
eff.plot(eff_Shannon, plotdata = T,
         ylab = "Number of effective pollen species on one pollinator",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_Shannon))
hist(resid(m_Shannon)) # residuals look good 

ks.test(resid(m_Shannon), "pnorm", mean = mean(resid(m_Shannon)), sd = sd(resid(m_Shannon)))
#significance, probably because of "tail" of residuals at higher values

residuals_Shannon <- simulateResiduals(fittedModel = m_Shannon)
plot(residuals_Shannon)
testOutliers(residuals_Shannon)
#KS test and combined adjusted quantile test significant, outliers non-significant


###quadratic----
#model ranked best with quadratic term (m7)
m_Shannon2 <- glmmTMB(Shannon_diversity ~ log(Stems) + I(log(Stems)^2) + Group 
                      + (1|Site) + (1|Species), data = dist_data)
summary(m_Shannon2)

eff_Shannon2 <- effect(c("log(Stems)","I(log(Stems)^2)"),m_Shannon2, xlevels = 50)  
eff.plot(eff_Shannon2, plotdata = T,
         ylab = "Number of effective pollen species on one pollinator",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)
#quadratic term is plotted, only visible on larger scale!

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_Shannon2))
hist(resid(m_Shannon2)) # residuals look okay

ks.test(resid(m_Shannon2), "pnorm", mean = mean(resid(m_Shannon)), sd = sd(resid(m_Shannon)))
#non-significant

residuals_Shannon2 <- simulateResiduals(fittedModel = m_Shannon2)
plot(residuals_Shannon2)
testOutliers(residuals_Shannon2)
#KS test and combined adjusted quantile test significant, outliers non-significant


#effect sizes----
r.squaredGLMM(m_Shannon)
#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems for pollinators caught on Arnica flowers and in surrouding area
pred_data_Shannon2 <- data.frame(Stems = rep(c(10, 100, 500),2), 
                                Group = c(rep("area",3),rep("flower",3)))
predictions_Shannon2 <- predict(m_Shannon2, newdata = pred_data_Shannon, type = "response", re.form = NA, se.fit = T)

pred_data_Shannon2$Pred.Poll_diver <- predictions_Shannon2$fit
pred_data_Shannon2$Pred.POll_diver_SE <- predictions_Shannon2$se.fit
pred_data_Shannon2
#In area:
#for 10 stems, 3.5 (+/- 0.2) effective species of pollen are expected to be carried
#by pollinators, for 100 stems 3.1 (+/- 0.1) effective species, and for 500 stems
#3.1 (+/- 0.2) effective species of pollen are expected to be carried by a pollinator.

#On flower:
#for 10 stems, 3.2 (+/- 0.2) effective species of pollen are expected to be carried
#by pollinators, for 100 stems 2.7 (+/- 0.1) effective species, and for 500 stems
#2.7 (+/- 0.2) effective species of pollen are expected to be carried by a pollinator.