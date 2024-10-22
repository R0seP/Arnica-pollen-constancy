#This script extracts pollen samples that consist of a maximum of 3 effective
#pollen species AND have more than 50% Arnica pollen. These samples are considered
#pollen constant on Arnica. The number of these constant samples is then modeled.

setwd("C:/Users/sohe1/documents/Master General Biology/Master_Thesis/R")
library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
library(vegan)
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/EffPlots.R")

#get data
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)

#prepare data----
adj <- read.csv("adj_new.csv", h = T)
pollen <- adj[,c(3,8:42)]
pollen <- na.omit(pollen)

#calculate Shannon (Renyi?) diversity of the pollen samples
distances <- renyi(pollen[,2:36], scales = c(0,1)) #calculate species richness and Shannon entropy
distances$filn <- pollen$filn
distances <- distances %>% select(,c("filn","0","1"))
distances <- distances %>%
  rename(Species_richness = "0", Shannon_entropy = "1")
distances$Shannon_diversity = 2^distances$Shannon_entropy #calculate Shannon diversity (as 2^H (H = Shannon entropy))
dist_data <- inner_join(distances, comb_all, by = "filn")

#add "Arnica constancy" 
dist_data$Arnica_constant <- rep(0,length(dist_data$filn))

for (i in 1:length(dist_data$filn)){
  dist_data[i, "Arnica_constant"] <- ifelse(dist_data[i,"Shannon_diversity"] <= 3 
                                      && dist_data[i, "P_ASTE.Arnica.montana"] > 0.5, 
                                      "yes", "no")
}     #Add Arnica constancy
head(dist_data$Arnica_constant)

Arnica_constant <- subset(dist_data, Arnica_constant == "yes") #subset only Arnica constant samples
Not_Arnica_constant <- subset(dist_data, Arnica_constant == "no") #subset only not Arnica constant samples

Arnica_constant_counts <- Arnica_constant %>%
  count(Site) #count the number of Arnica constant species per site
Not_Arnica_constant_counts <- Not_Arnica_constant %>%
  count(Site) #count the nuber of not Arnica constant samples per site

Arnica_constant_counts <- Arnica_constant_counts %>%
  rename(n_Arnica_constant = n) #rename column with counts for Arnica constant samples
Not_Arnica_constant_counts <- Not_Arnica_constant_counts %>%
  rename(n_Not_Arnica_constant = n) #rename column with counts for not Arnica constant samples

constancy <- inner_join(Arnica_constant_counts, Not_Arnica_constant_counts, 
                        by = "Site") #combine counts for Arnica constant and not Arnica constant
constancy$Stems <- comb1$Stems #add population size


#site model selection----
m1 <- glmmTMB(cbind(n_Arnica_constant, n_Not_Arnica_constant) ~ Stems,
                      data = constancy, family = binomial)
m2 <- glmmTMB(cbind(n_Arnica_constant, n_Not_Arnica_constant) ~ log(Stems),
              data = constancy, family = binomial)

mlist = list(m1, m2)
AICTab = AIC(m1, m2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
#model m2 with log ranked higher

r.squaredGLMM(m1)
r.squaredGLMM(m2)
#model m2 with log also higher r2, use that model!

#site model----
#Is the number of samples constant on Arnica influenced by the population size?

m_constancy <- glmmTMB(cbind(n_Arnica_constant, n_Not_Arnica_constant) ~ log(Stems),
                       data = constancy, family = binomial)
summary(m_constancy)
r.squaredGLMM(m_constancy)

check_overdispersion(m_constancy)
#no overdispersion

eff_constancy <- effect("log(Stems)",m_constancy, xlevels = 50)  
eff.plot(eff_constancy, plotdata = T,
         ylab = "Proportion of Arnica-constant samples",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_constancy))
hist(resid(m_constancy))
#qqplot looks good, hist of residuals looks more uniformly distributed

ks.test(resid(m_constancy), "pnorm", mean = mean(resid(m_constancy)), sd = sd(resid(m_constancy)))
#ks test not significant

residuals_constancy <- simulateResiduals(fittedModel = m_constancy)
plot(residuals_constancy)
testOutliers(residuals_constancy)
#no significant deviations or outliers

#effect sizes site model----

#r-squared value:
r.squaredGLMM(m_constancy)

#predict proportion of filled seeds for 10, 100, and 500 stems
pred_data_constancy <- data.frame(Stems = c(10, 100, 500))
predictions_constancy <- predict(m_constancy, newdata = pred_data_constancy, type = "response", re.form = NA, se.fit = T)

pred_data_constancy$Pred.Arnica.constant <- predictions_constancy$fit
pred_data_constancy$Pred.Arnica.constant_SE <- predictions_constancy$se.fit
pred_data_constancy
#for 10 stems, 16.0 (+/- 2.8)% of sample are expected to be Arnica-constant, for 
#100 stems 24.4 (+/- 1.9)% of samples, and for 500 stems 31.8 (+/- 3.8)% of 
#samples are expected to be constant on Arnica (have more than 50% Arnica and
#up to three effective pollen species.


#prepare data 2----
dist_data$Arnica_constant_num <- rep(0,length(dist_data$filn))
for (i in 1:length(dist_data$filn)){
  dist_data[i, "Arnica_constant_num"] <- ifelse(dist_data[i,"Arnica_constant"] == "yes", 
                                                1, 0)
}     #Add column with 1 where Arnica constant and 0 where not Arnica constant
head(dist_data$Arnica_constant_num)

constancy2 <- na.omit(dist_data)


#all samples model selection----
m1 <- glmmTMB(Arnica_constant_num ~ Stems + (1|Species) + (1|Site),
              family = binomial, data = constancy2)
m2 <- glmmTMB(Arnica_constant_num ~ log(Stems) + (1|Species) + (1|Site),
              family = binomial, data = constancy2)

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
#r2 as NA at the moment

#all samples model----
m_constancy2 <- glmmTMB(Arnica_constant_num ~ Stems + (1|Species) + (1|Site),
              family = binomial, data = constancy2)
summary(m_constancy2)

eff_constancy2 <- effect("Stems",m_constancy2, xlevels = 50)  
eff.plot(eff_constancy2, plotdata = T,
         ylab = "Proportion of Arnica-constant samples",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)

#test if model assumptions are met and test model for fit:
qqnorm(resid(m_constancy2))
hist(resid(m_constancy2))
#qqplot and histogram look horrible, maybe go with first model
check_overdispersion(m_constancy2)
#no overdispersion

resid_jittered <- resid(m_constancy2) + rnorm(length(resid(m_constancy2)), mean = 0, sd = 1e-6)
ks.test(resid_jittered, "pnorm", mean = mean(resid_jittered), sd = sd(resid_jittered))
#ks test significant

residuals_constancy2 <- simulateResiduals(fittedModel = m_constancy2)
plot(residuals_constancy2)
testOutliers(residuals_constancy2)
#no significant deviations or outliers, weird because KS test significant before?

#effect sizes all samples----
#predict proportion of filled seeds for 10, 100, and 500 stems
pred_data_constancy2 <- data.frame(Stems = c(10, 100, 500))
predictions_constancy2 <- predict(m_constancy2, newdata = pred_data_constancy2, type = "response", re.form = NA, se.fit = T)

pred_data_constancy2$Pred.Arnica.constant <- predictions_constancy2$fit
pred_data_constancy2$Pred.Arnica.constant_SE <- predictions_constancy2$se.fit
pred_data_constancy2
#for 10 stems, 10.5 (+/- 3.5)% of sample are expected to be Arnica-constant, for 
#100 stems 11.5 (+/- 3.7)% of samples, and for 500 stems 16.9 (+/- 6.5)% of 
#samples are expected to be constant on Arnica (have more than 50% Arnica and
#up to three effective pollen species.

#effect sizes from second model quite different from first model
#first model m_constancy seems to have better fit
