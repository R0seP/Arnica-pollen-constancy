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
min(comb_all2$P_ASTE.Arnica.montana)  #0%
max(comb_all2$P_ASTE.Arnica.montana)  #98.1%

#Nr of Arnica pollen carried 
mean(comb_all2$Nr_Arnica) #273
min(comb_all2$Nr_Arnica)  #0
max(comb_all2$Nr_Arnica)  #2554
