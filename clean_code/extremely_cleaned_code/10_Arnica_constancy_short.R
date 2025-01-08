#This script extracts pollen samples that consist of a maximum of 3 effective
#pollen species AND have more than 50% Arnica pollen. These samples are considered
#pollen constant on Arnica. The number of these constant samples is then modeled.

library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(effects)
library(DHARMa)
library(MuMIn)
library(performance)
library(vegan)
source("EffPlots.R") #for the code for "EffPlots.R" please contact Ola Olsson, Lund University, Biology

#get data
adj <- read.csv("adj_new.csv", h = T)
comb_all <- read.csv("comb_all.csv", h = T)
comb1 <- read.csv("comb1.csv", h = T)

#prepare data----
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
}    #Add Arnica constancy
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


#model----
#Is the number of samples constant on Arnica influenced by the population size?

m_constancy <- glmmTMB(cbind(n_Arnica_constant, n_Not_Arnica_constant) ~ log(Stems),
                       data = constancy, family = binomial)
summary(m_constancy)

#r-squared value:
r.squaredGLMM(m_constancy)

eff_constancy <- effect("log(Stems)",m_constancy, xlevels = 50)  
eff.plot(eff_constancy, plotdata = T,
         ylab = "Proportion of Arnica-constant samples",
         xlab = "Population size Arnica (Nr Stems)",
         main = "",
         ylim.data = T, overlay = F, col.data = 3)


#test if model assumptions are met and test model for fit:
check_overdispersion(m_constancy)
#no overdispersion
hist(resid(m_constancy))

ks.test(resid(m_constancy), "pnorm", mean = mean(resid(m_constancy)), sd = sd(resid(m_constancy)))
#ks test not significant

residuals_constancy <- simulateResiduals(fittedModel = m_constancy)
plot(residuals_constancy)
testOutliers(residuals_constancy)
#no significant deviations or outliers

#effect sizes----
#predict proportion of constant samples for 10, 100, and 500 stems
pred_data_constancy <- data.frame(Stems = c(10, 100, 500))
predictions_constancy <- predict(m_constancy, newdata = pred_data_constancy, type = "response", re.form = NA, se.fit = T)

pred_data_constancy$Pred.Arnica.constant <- predictions_constancy$fit
pred_data_constancy$Pred.Arnica.constant_SE <- predictions_constancy$se.fit
pred_data_constancy
#for 10 stems, 16.0 (+/- 2.8)% of samples are expected to be Arnica-constant, for 
#100 stems 24.4 (+/- 1.9)% of samples, and for 500 stems 31.8 (+/- 3.8)% of 
#samples are expected to be constant on Arnica (have more than 50% Arnica and
#up to three effective pollen species.