#This script summarizes all models used in the thesis.

library(glmmTMB)
library(lmerTest)


#model direct influence of population size on reproductive success
seed_data <- read.csv("seed_data.csv", h = T)

##proportion of filled seeds as response
m_success2 <- glmmTMB(Filled ~ log(Stems) + offset(log(Total_seeds))  + (1|Site),
                      data = seed_data, family = nbinom2)

##number of filled seeds as response
m_success4 <- glmmTMB(Filled ~ log(Stems)  + (1|Site),
                      data = seed_data, family = nbinom1)


#model number of pollinators
poll_nr <- read.csv("poll_nr.csv", h = T)

poll_nr$Date_from_start_factor <- as.factor(poll_nr$Date_from_start)
poll_nr$Poll_per_hr <- round(poll_nr$Poll_per_hr)

m_poll_nr <- glmmTMB(Poll_per_hr ~ Stems + Date_from_start
                     + (1|time_quali) + (1|Site), data = poll_nr, family = poisson)


#model proportion of A. montana associated pollinators
comb1 <- read.csv("comb1.csv", h = T)

m_community2 <- glmmTMB(cbind(Arnica_associated_noMel, Area_associated) ~ log(Stems), 
                        data = comb1, family = binomial)


#model proportion of A. montana pollen carried
comb_all2 <- read.csv("comb_all2.csv", h = T)

m_species <- glmmTMB(Nr_Arnica ~ Species + Group * Stems + offset(log(nPoll)) + (1|Site), 
                     data = comb_all2, family = nbinom1,
                     control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000)))


#model specialization index Blüthgen's d
species_metrics <- read.csv("species_metrics.csv", h = T)

m_d <-  glmmTMB(Arnica_delta_d ~ log(Stems), data = species_metrics)


#model Nr. of A. montana constant samples
constancy <- read.csv("constancy.csv", h = T)

m_constancy <- glmmTMB(cbind(n_Arnica_constant, n_Not_Arnica_constant) ~ log(Stems),
                       data = constancy, family = binomial)


#model Nr. of area-caught pollinators that carried >= 5% A. montana pollen
Arnica5_counts <- read.csv("Arnica5_counts.csv", h = T)

m_Arnica5_area <- glmmTMB(n ~ Stems, data = Arnica5_counts, family = poisson)


