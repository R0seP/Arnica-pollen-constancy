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
adj <- read.csv("pollen_adj_class.csv", h = T)


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
m5 <- glmmTMB(Shannon_diversity ~ log(Stems) + (1|Site) + (1|Species), 
              data = dist_data, family = poisson)

mlist = list(m1, m2, m3, m4, m5)
AICTab = AIC(m1, m2, m3, m4, m5) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#m4 without Group, with gaussian family and with log ranked highest

r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
r.squaredGLMM(m4)
#r2 of m2 highest, still use without interaction to use simple model 

#model----
#Model the effect of Arnica population size on the effective number of species a 
#pollinator was carrying

m_Shannon <- glmmTMB(Shannon_diversity ~ log(Stems) + (1|Site) + (1|Species), 
                     data = dist_data)
summary(m_Shannon)

m_Shannon2 <- glmmTMB(log(Shannon_diversity) ~ log(Stems) + I(log(Stems)^2) + (1|Site) + (1|Species), 
                     data = dist_data)
summary(m_Shannon2)

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

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_Shannon2))
hist(resid(m_Shannon2)) # residuals look good 

ks.test(resid(m_Shannon2), "pnorm", mean = mean(resid(m_Shannon)), sd = sd(resid(m_Shannon)))
#significance, probably because of "tail" of residuals at higher values

residuals_Shannon2 <- simulateResiduals(fittedModel = m_Shannon2)
plot(residuals_Shannon2)
testOutliers(residuals_Shannon2)
#KS test and combined adjusted quantile test significant, outliers non-significant


#effect sizes----
r.squaredGLMM(m_Shannon)
#predict proportion of Arnica-associated pollinators in community
#for 10, 100, and 500 stems
pred_data_Shannon <- data.frame(Stems = c(10, 100, 500))
predictions_Shannon <- predict(m_Shannon, newdata = pred_data_Shannon, type = "response", re.form = NA, se.fit = T)

pred_data_Shannon$Pred.Poll_diver <- predictions_Shannon$fit
pred_data_Shannon$Pred.POll_diver_SE <- predictions_Shannon$se.fit
pred_data_Shannon
#for 10 stems, 3.3 (+/- 0.1) effective species of pollen are expected to be carried
#by pollinators, for 100 stems 3.0 (+/- 0.1) effective species, and for 500 stems
#2.8 (+/- 0.1) effective species of pollen are expected to be carried by a pollinator.