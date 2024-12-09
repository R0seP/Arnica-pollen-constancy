##this is a script for exploring my complete data set and trying things out
##in the analysis

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")
library(tidyverse)
library(glmmTMB)
library(lmerTest)
library(gridExtra)
library(effects)
library(lattice)
library(permute)
library(vegan)

#get data
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)


######################################################################################################
#explore relationship reproductive success ~ population size
plot(Av_filling_rate ~ Stems, data = comb1)
#looks as though there could be a trend for more filling rate the more stems
#relationship looks logarithmic, transform?
plot(Av_filling_rate ~ log(Stems), data = comb1)

summary(lm(Av_filling_rate ~ log(Stems), data = comb1))
#no significance, positive trend
summary(lm(Av_filling_rate ~ Stems, data = comb1))

ggplot(comb1, aes(x = log(Stems), y = Av_filling_rate)) +
  geom_point() +  # Add data points
  geom_smooth(method = "lm", se = T)
#cannot isolate genetic component since influence of population size will always include 
#both genetic and environmental (e.g., pollinator) effects?


###########################################################################################################
#Explore data with further graphs----

boxplot(comb_all$P_ASTE.Arnica.montana~comb_all$Group)
#much more Arnica in flower group than in area group

boxplot(comb_all$P_ASTE.Arnica.montana~comb_all$Overview)
#a lot of Arnica pollen on Beetle and Fly, some on Moth, Butterfly and Bee, 
#barely any on Sawfly and Wasp
boxplot(comb_all$P_ASTE.Arnica.montana~comb_all$Family, par(las = 2))
boxplot(comb_all$P_ASTE.Arnica.montana~comb_all$Species, par(las = 2))
#some species (Eristalis, Merodon, Empis, Dasytes, ...) much more Arnica than others

plot(comb_all$P_ASTE.Arnica.montana~comb_all$Stems)
plot(comb_all[comb_all$Group=="flower", ]$P_ASTE.Arnica.montana~comb_all[comb_all$Group=="flower", ]$Stems)
plot(comb_all[comb_all$Group=="area", ]$P_ASTE.Arnica.montana~comb_all[comb_all$Group=="area", ]$Stems)
#very uniform distribution of Arnica carried in both groups with the population size

plot(comb_all$Av_filling_rate~comb_all$P_ASTE.Arnica.montana)
#very uniformly distributed

plot(comb1$Eristalis~comb1$Stems)
#pretty scattered, maybe increasing Eristalis with more stems?

plot(comb1$Av_filling_rate~comb1$Eristalis)
#tendentially increasing filling rate with more Eristalis

ggplot(comb_all, aes(fill=Group, y=Nr_Arnica, x=Group)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Sampling group", y = "Nr of Arnica montana pollen carried")

ggplot(comb_all, aes(fill=Group, y=P_ASTE.Arnica.montana, x=Group)) + 
  geom_violin()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Sampling group", y = "% Arnica montana pollen carried")

ggplot(comb_all2, aes(fill=Group, y=P_ASTE.Arnica.montana, x=Species)) + 
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
  scale_fill_manual(values = c("area" = "darkgreen","flower"="orange"))+
  labs(x = "Species", y = "% Arnica montana pollen carried")

mosaicplot(Species ~ Group, data = comb_all2,
           las = 2,
           main = "",
           ylab = "Group",
           xlab = "Pollinator species",
           col = c("darkgreen","orange"),
           cex = 0.9)

plot(Av_filling_rate ~ Stems, data = comb1)
plot(plogis(Av_filling_rate) ~ Stems, data = comb1) #plogis(0.5) calculates the 
#probability that a logistic random variable is less than or equal to 0.5.
####################################################################################################
#what influences the %Arnica carried?###

#m_pollen <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Species + Stems +
#                    (1|Site) + (1+Species|Group), family = binomial, 
#                  na.action = na.omit, data = comb_all)

#summary(m_pollen)
#qqnorm(resid(m_pollen))
#model leads to fatal error in R. Probably not ideal anyways since random factor
#"Group" only has two levels



#predictions: 
comb_all2$predictions <- predict(m_pollen, type = "response")

ggplot(comb_all2, aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "pred. m_pollen", x = "  Nr Stems", y = "Predicted % Arnica carried", color = "Group")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

#by species:
p1 <- ggplot(comb_all2[comb_all2$Species == "Andrena sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Andrena sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p2 <- ggplot(comb_all2[comb_all2$Species == "Apis mellifera",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Apis mellifera", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p3 <- ggplot(comb_all2[comb_all2$Species == "Bombus pascuorum",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Bombus pascuorum", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p4 <- ggplot(comb_all2[comb_all2$Species == "Bombus ruderarius",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Bombus ruderarius", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p5 <- ggplot(comb_all2[comb_all2$Species == "Bombus terrestris",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Bombus terrestris", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p6 <- ggplot(comb_all2[comb_all2$Species == "Coenonympha pamphilus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Coenonympha pamphilus", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p7 <- ggplot(comb_all2[comb_all2$Species == "Dasytes niger",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Dasytes niger", x = "", y = "Predicted % Arnica carried", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p8 <- ggplot(comb_all2[comb_all2$Species == "Empis livida",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Empis livida", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p9 <- ggplot(comb_all2[comb_all2$Species == "Empis tessellata",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Empis tesellata", x = " ", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p10 <- ggplot(comb_all2[comb_all2$Species == "Eristalis sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Eristalis sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p11 <- ggplot(comb_all2[comb_all2$Species == "Eupeodes corollae",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Eupeodes corollae", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p12 <- ggplot(comb_all2[comb_all2$Species == "Helophilus pendulus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Helophilus pendulus", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p13 <- ggplot(comb_all2[comb_all2$Species == "Lasioglossum sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Lasioglossum sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p14 <- ggplot(comb_all2[comb_all2$Species == "Maniola jurtina",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Maniola jurtina", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p15 <- ggplot(comb_all2[comb_all2$Species == "Meligethes sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Meligethes sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p16 <- ggplot(comb_all2[comb_all2$Species == "Merodon equestris",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Merodon equestris", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p17 <- ggplot(comb_all2[comb_all2$Species == "Nomada sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Nomada sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p18 <- ggplot(comb_all2[comb_all2$Species == "Ochlodes sylvanus",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Ochlodes sylvanus", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p19 <- ggplot(comb_all2[comb_all2$Species == "Oedemera sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Oedemera sp", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p20 <- ggplot(comb_all2[comb_all2$Species == "Phyllopertha horticola",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Phyllopertha horticola", x = "", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p21 <- ggplot(comb_all2[comb_all2$Species == "Sphaerophoria sp",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Sphaerophoria sp", x = "Nr Stems", y = "", color = "Group")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

p22 <- ggplot(comb_all2[comb_all2$Species == "Stenurella melanura",], aes(x = Stems, y = predictions, color = as.factor(Group))) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Stenurella melanura", x = "", y = "", color = "Group") +
  theme(legend.position = c(1.2, 0.5))+
  scale_color_manual(values = c("area" = "darkgreen","flower"="orange"))

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
             p15, p16, p17, p18, p19, p20, p21, p22, ncol = 6)


##############################################################################################################
#is the number of frequent species influenced by the number of stems?

#important (according to simper) species that prefer Arnica
plot(Eristalis_sp ~ Stems, data = comb1)
n_Eristalis <- glmmTMB(Eristalis_sp ~ Stems, family = poisson, data = comb1)
summary(n_Eristalis)

plot(Meligethes_sp ~ Stems, data = comb1)
n_Meligethes <- glmmTMB(Meligethes_sp ~ Stems, family = poisson, data = comb1)
summary(n_Meligethes)

plot(Stenurella_melanura ~ Stems, data = comb1)
n_Stenurella <- glmmTMB(Stenurella_melanura ~ Stems, data = comb1, family = poisson)
summary(n_Stenurella)

plot(Dasytes_niger ~ Stems, data = comb1)
n_Dasytes <- glmmTMB(Dasytes_niger ~ Stems, data = comb1, family = poisson)
summary(n_Dasytes)

plot(Empis_livida ~ Stems, data = comb1)
n_Empis_livida <- glmmTMB(Empis_livida ~ Stems, data = comb1, family = poisson)
summary(n_Empis_livida)

plot(Helophilus_pendulus ~ Stems, data = comb1)
n_Helophilus <- glmmTMB(Helophilus_pendulus ~ Stems, data = comb1, family = poisson)
summary(n_Helophilus)

plot(Empis_tessellata ~ Stems, data = comb1)
n_Empis_tessellata <- glmmTMB(Empis_tessellata ~ Stems, data = comb1, family = poisson)
summary(n_Empis_tessellata)

plot(Oedemera_sp ~ Stems, data = comb1)
n_Oedemera <- glmmTMB(Oedemera_sp ~ Stems, data = comb1, family = poisson)
summary(n_Oedemera)

#important (according to simper) species that prefer flowers other than Arnica:
plot(Bombus_terrestris ~ Stems, data = comb1)
n_Bombus_terrestris <- glmmTMB(Bombus_terrestris ~ Stems, data = comb1, family = poisson)
summary(n_Bombus_terrestris)

plot(Bombus_ruderarius ~ Stems, data = comb1)
n_Bombus_ruderarius <- glmmTMB(Bombus_ruderarius ~ Stems, data = comb1, family = poisson)
summary(n_Bombus_ruderarius)

plot(Apis_mellifera ~ Stems, data = comb1)
n_Apis <- glmmTMB(Apis_mellifera ~ Stems, data = comb1, family = poisson)
summary(n_Apis)

#Syrphid species that seems very generalized (no group preference according to simper):
plot(Eupeodes_corollae ~ Stems, data = comb1)
n_Eupeodes <- glmmTMB(Eupeodes_corollae ~ Stems, data = comb1, family = poisson)
summary(n_Eupeodes)

#additional Syrphid that seems important for Arnica pollination:
plot(Merodon_equestris ~ Stems, data = comb1)
n_Merodon <- glmmTMB(Merodon_equestris ~ Stems, data = comb1, family = poisson)
summary(n_Merodon)

#################################################################################################################
#is the amount of pollen carried in the important species dependent on the population size of Arnica?

subset_Eristalis <- subset(comb_all2, Species == "Eristalis sp")
m_Eristalis <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                      (1|Site), family = binomial, 
                    na.action = na.omit, data = subset_Eristalis)
summary(m_Eristalis)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Eristalis, main = "Eristalis, all samples")


subset_Meligethes <- subset(comb_all2, Species == "Meligethes sp")
m_Meligethes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Meligethes)
summary(m_Meligethes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Meligethes, main = "Meligethes, all samples")


subset_Stenurella <- subset(comb_all2, Species == "Stenurella melanura")
m_Stenurella <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Stenurella)
summary(m_Stenurella)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Stenurella, main = "Stenurella, all samples")


subset_Dasytes <- subset(comb_all2, Species == "Dasytes niger")
m_Dasytes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Dasytes)
summary(m_Dasytes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Dasytes, main = "Dasytes, all samples")


subset_Empis_livida <- subset(comb_all2, Species == "Empis livida")
m_Empis_livida <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Empis_livida)
summary(m_Empis_livida)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Empis_livida, main = "Empis livida, all samples")


subset_Empis_tessellata <- subset(comb_all2, Species == "Empis tessellata")
m_Empis_tessellata <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Empis_tessellata)
summary(m_Empis_tessellata)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Empis_tessellata, main = "Empis tessellata, all samples")


subset_Helophilus <- subset(comb_all2, Species == "Helophilus pendulus")
m_Helophilus <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Helophilus)
summary(m_Helophilus)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Helophilus, main = "Helophilus, all samples")


subset_Oedemera <- subset(comb_all2, Species == "Oedemera sp")
m_Oedemera <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Oedemera)
summary(m_Oedemera)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Oedemera, main = "Oedemera, all samples")


subset_Bombus_terrestris <- subset(comb_all2, Species == "Bombus terrestris")
m_Bombus_terrestris <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Bombus_terrestris)
summary(m_Bombus_terrestris)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Bombus_terrestris, main = "Bombus terrestris, all samples")


subset_Bombus_ruderarius <- subset(comb_all2, Species == "Bombus ruderarius")
m_Bombus_ruderarius <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Bombus_ruderarius)
summary(m_Bombus_ruderarius)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Bombus_ruderarius, main = "Bombus ruderarius, all samples")


subset_Apis <- subset(comb_all2, Species == "Apis mellifera")
m_Apis <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Apis)
summary(m_Apis)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Apis, main = "Apis, all samples")


subset_Merodon <- subset(comb_all2, Species == "Merodon equestris")
m_Merodon <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Merodon)
summary(m_Merodon)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Merodon, main = "Merodon, all samples")


subset_Eupeodes <- subset(comb_all2, Species == "Eupeodes corollae")
m_Eupeodes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems * Group +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Eupeodes)
summary(m_Eupeodes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Eupeodes, main = "Eupeodes, all samples")


#########################################################################################################
#look at amount of pollen carried only in the flower-insects

subset_Eristalis <- subset(subset_Eristalis, Group == "flower")
m_Eristalis <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                         (1|Site), family = binomial, 
                       na.action = na.omit, data = subset_Eristalis)
summary(m_Eristalis)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Eristalis, main = "Eristalis, flower group")


subset_Meligethes <- subset(subset_Meligethes, Group == "flower")
m_Meligethes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                          (1|Site), family = binomial, 
                        na.action = na.omit, data = subset_Meligethes)
summary(m_Meligethes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Meligethes, main = "Meligethes, flower group")


subset_Stenurella <- subset(subset_Stenurella, Group == "flower")
m_Stenurella <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                          (1|Site), family = binomial, 
                        na.action = na.omit, data = subset_Stenurella)
summary(m_Stenurella)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Stenurella, main = "Stenurella, flower group")


subset_Dasytes <- subset(subset_Dasytes, Group == "flower")
m_Dasytes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                       (1|Site), family = binomial, 
                     na.action = na.omit, data = subset_Dasytes)
summary(m_Dasytes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Dasytes, main = "Dasytes, flower group")


subset_Empis_livida <- subset(subset_Empis_livida, Group == "flower")
m_Empis_livida <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                            (1|Site), family = binomial, 
                          na.action = na.omit, data = subset_Empis_livida)
summary(m_Empis_livida)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Empis_livida, main = "Empis livida, flower group")


subset_Empis_tessellata <- subset(subset_Empis_tessellata, Group == "flower")
m_Empis_tessellata <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                (1|Site), family = binomial, 
                              na.action = na.omit, data = subset_Empis_tessellata)
summary(m_Empis_tessellata)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Empis_tessellata, main = "Empis tessellata, flower group")


subset_Helophilus <- subset(subset_Helophilus, Group == "flower")
m_Helophilus <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                          (1|Site), family = binomial, 
                        na.action = na.omit, data = subset_Helophilus)
summary(m_Helophilus)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Helophilus, main = "Helophilus, flower group")


subset_Oedemera <- subset(subset_Oedemera, Group == "flower")
m_Oedemera <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                        (1|Site), family = binomial, 
                      na.action = na.omit, data = subset_Oedemera)
summary(m_Oedemera)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Oedemera, main = "Oedemera, flower group")


subset_Bombus_terrestris <- subset(subset_Bombus_terrestris, Group == "flower")
m_Bombus_terrestris <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                                 (1|Site), family = binomial, 
                               na.action = na.omit, data = subset_Bombus_terrestris)
summary(m_Bombus_terrestris)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Bombus_terrestris, main = "Bombus terrestris, flower group")


#subset_Bombus_ruderarius <- subset(subset_Bombus_ruderarius, Group == "flower")
#m_Bombus_ruderarius <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
#                                 (1|Site), family = binomial, 
#                               na.action = na.omit, data = subset_Bombus_ruderarius)
#summary(m_Bombus_ruderarius)
#plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Bombus_ruderarius, , main = "Bombus ruderarius, flower group")
#Bombus ruderarius only caught in area!


subset_Apis <- subset(subset_Apis, Group == "flower")
m_Apis <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                    (1|Site), family = binomial, 
                  na.action = na.omit, data = subset_Apis)
summary(m_Apis)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Apis, main = "Apis, flower group")


subset_Merodon <- subset(subset_Merodon, Group == "flower")
m_Merodon <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                       (1|Site), family = binomial, 
                     na.action = na.omit, data = subset_Merodon)
summary(m_Merodon)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Merodon, main = "Merodon, flower group")


subset_Eupeodes <- subset(subset_Eupeodes, Group == "flower")
m_Eupeodes <- glmmTMB(cbind(Nr_Arnica, Nr_Not.Arnica) ~ Stems +
                        (1|Site), family = binomial, 
                      na.action = na.omit, data = subset_Eupeodes)
summary(m_Eupeodes)
plot(P_ASTE.Arnica.montana ~ Stems, data = subset_Eupeodes, main = "Eupeodes, flower group")



#######################################################################################################
#Model the reproductive success of Arnica

m_success <- glmmTMB(Av_filling_rate ~ log(Stems) + Species + P_ASTE.Arnica.montana +
                       (1|Site), data = comb_all)
summary(m_success)
qqnorm(resid(m_success)) #residuals look terrible without a log, maybe fine with log?

###very confusing, negative impact of more Arnica pollen & bigger population size on reproductive success
