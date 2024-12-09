#this script calculates overview statistics over the pollen grains carried in 
#the Arnica study

library(tidyverse)
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)
comb_all2 <- na.omit(comb_all2)

#number of pollen carried
mean(comb_all2$nPoll) #1095
min(comb_all2$nPoll)  #1
max(comb_all2$nPoll)  #21223

#% of Arnica carried
mean(comb_all2$P_ASTE.Arnica.montana) #28.4%
mean(comb_all2[comb_all2$Group == "flower","P_ASTE.Arnica.montana"]) # 49.4%
mean(comb_all2[comb_all2$Group == "area","P_ASTE.Arnica.montana"]) # 4.2%

min(comb_all2$P_ASTE.Arnica.montana)  #0%
max(comb_all2$P_ASTE.Arnica.montana)  #98.1%

median(comb_all2$P_ASTE.Arnica.montana) #10.2%
median(comb_all2[comb_all2$Group == "flower","P_ASTE.Arnica.montana"]) # 50.1%
median(comb_all2[comb_all2$Group == "area","P_ASTE.Arnica.montana"]) # 0.4%

#Nr of Arnica pollen carried 
mean(comb_all2$Nr_Arnica) #273
min(comb_all2$Nr_Arnica)  #0
max(comb_all2$Nr_Arnica)  #2554

#Nr of samples with more than 95% Arnica
Arnica_constant_95 <- comb_all2[which(comb_all2$P_ASTE.Arnica.montana >= 0.95), ]
