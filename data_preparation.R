#This script takes the various tables containing field data and lab analysis results
#and combines everything into several data frames.

#The data frame "comb_all" (n = 568) consists of the pollen identification results
#with adjusted classification for every sample, while also including population
#size and average seed filling rate at the site the pollinator was caught at. It
#also includes the species and the group (caught on Arnica flower or in surrounding
#area) for each pollinator.

#The data frame "comb_all2" (n = 449 is a subset of "comb_all" and includes the same
#data as comb_all, but only for species that were observed at least 5 times (found
#in "imp_species").

#The data frame "comb1" (n = 20) shows the population size and the average reproductive
#success of Arnica at each site, it also includes counts of species that were
#caught five times or more (to be found in imp_species as well), and the categories
#"Arnica_associated" and "Area_associated", which adds up all species that have 
#been found to favor Arncia or flowers other than Arnica respectively. 

#The data frames "Arnica_polli" (n = 226) and "Area_polli" (n = 218) are two subsets  
#of "comb_all2" that include only samples of species that have been shown to be associated
#more with Arnica or the surrounding area respectively.

#The data frame "seed_data" (n = 381) contains the seed counts for every flower
#head combined with population size.

#The data frame "poll_nr" (n = 41) combines the number of pollinators caught at each 
#sampling occasion with the population size of Arnica and meteorological as well
#as temporal data. 

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")

##import data
ori <- read.csv("pollen_orig_class.csv", h = T) #original classifications for pollen identity
adj <- read.csv("adj_new.csv", h = T)  #adjusted classifications for pollen identity
polli <- read.csv("pollinators.csv", h = T)     #information about pollinator species and group of samples
pop <- read.csv("pop_size.csv", h = T)          #population size of Arnica at the sites
seeds <- read.csv("seed_counts.csv", h = T)     #seed counts at sites as proxy for reproductive success


#get average filling ratio for all sites
av_seeds <- data.frame(Site = 1:20, Av_filling_rate = rep(0,20))

for (i in 1:20) {
  av_seeds[av_seeds$Site == i, "Av_filling_rate"] <- mean(seeds[seeds$Site == i, "Filling_Ratio"], 
                                                          na.rm = TRUE)}

#combine reproductive success and population size
library(dplyr)
comb1 <- inner_join(pop,av_seeds, by="Site")

#get bigger data frame with all seed counts and population size
seed_data <- inner_join(pop,seeds, by="Site")
seed_data$Total_seeds <- seed_data$Filled + seed_data$Not_Filled
seed_data$Site <- as.factor(seed_data$Site)

#combine pollinator and pollen identity data 
comb2 <- inner_join(adj, polli, by = c("Site", "Pollinator"))

#add population size data and reproductive success data:
comb_all <- comb2 %>%
  left_join(pop, by = "Site") %>%
  left_join(av_seeds, by = "Site")

comb_all <- comb_all[,-48] #delete comments column

comb_all <- comb_all[, c(1:3, 43:49, 4:42)] #re-sort for better readability

#create column with percentage of all other pollen
comb_all$P_Not.Arnica <- rep(0, length(comb_all$Site))
for (i in 1:length(comb_all$Site)){
  comb_all[i, "P_Not.Arnica"] <- sum(comb_all[i, c(15:16,18:49)])
}

library(this.path)
script_path <- this.path()
print(script_path)

#extract species with more than 5 occurrences
species_counts <- comb_all%>%group_by(Species)%>%summarise(Count = n())
imp_species <- species_counts[species_counts$Count >= 5, ]
imp_species <- imp_species[-23, ]

result <- comb_all %>%
  group_by(Site, Species) %>%
  summarise(Count = n())

Andrena_sp <- data.frame(Site = 1:20, Andrena_sp = rep(0,20))
result_Andrena<- na.omit(result[result$Species == "Andrena sp", ])
print(result_Andrena, n = 20)
Andrena_sp$Andrena_sp <- c(1,0,4,0,1,1,0,0,0,2,0,0,2,0,0,0,3,0,0,1)

Apis_mellifera <- data.frame(Site = 1:20, Apis_mellifera = rep(0,20))
result_Apis_mellifera<- na.omit(result[result$Species == "Apis mellifera", ])
print(result_Apis_mellifera, n = 20)
Apis_mellifera$Apis_mellifera <- c(2,0,1,0,4,2,7,0,5,1,1,4,15,0,2,0,2,4,0,0)

Bombus_pascuorum <- data.frame(Site = 1:20, Bombus_pascuorum = rep(0,20))
result_Bombus_pascuorum<- na.omit(result[result$Species == "Bombus pascuorum", ])
print(result_Bombus_pascuorum, n = 20)
Bombus_pascuorum$Bombus_pascuorum <- c(0,0,0,0,1,0,0,0,1,0,1,2,0,0,0,1,0,0,1,1)

Bombus_ruderarius <- data.frame(Site = 1:20, Bombus_ruderarius = rep(0,20))
result_Bombus_ruderarius<- na.omit(result[result$Species == "Bombus ruderarius", ])
print(result_Bombus_ruderarius, n = 20)
Bombus_ruderarius$Bombus_ruderarius <- c(0,1,0,0,0,0,1,0,4,1,0,0,0,2,0,0,0,1,0,0)

Bombus_terrestris <- data.frame(Site = 1:20, Bombus_terrestris = rep(0,20))
result_Bombus_terrestris<- na.omit(result[result$Species == "Bombus terrestris", ])
print(result_Bombus_terrestris, n = 20)
Bombus_terrestris$Bombus_terrestris <- c(0,1,0,0,0,0,1,0,1,0,5,1,1,0,1,0,3,2,0,1)

Coenonympha_pamphilus <- data.frame(Site = 1:20, Coenonympha_pamphilus = rep(0,20))
result_Coenonympha_pamphilus<- na.omit(result[result$Species == "Coenonympha pamphilus", ])
print(result_Coenonympha_pamphilus, n = 20)
Coenonympha_pamphilus$Coenonympha_pamphilus <- c(1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1)

Dasytes_niger <- data.frame(Site = 1:20, Dasytes_niger = rep(0,20))
result_Dasytes_niger<- na.omit(result[result$Species == "Dasytes niger", ])
print(result_Dasytes_niger, n = 20)
Dasytes_niger$Dasytes_niger <- c(3,0,2,4,2,0,2,0,3,0,0,1,0,0,0,0,1,0,6,0)

Empis_livida <- data.frame(Site = 1:20, Empis_livida = rep(0,20))
result_Empis_livida<- na.omit(result[result$Species == "Empis livida", ])
print(result_Empis_livida, n = 20)
Empis_livida$Empis_livida <- c(0,0,0,0,0,0,0,2,0,3,3,1,0,0,0,0,0,1,1,1)
                               
Empis_tessellata <- data.frame(Site = 1:20, Empis_tessellata = rep(0,20))
result_Empis_tessellata<- na.omit(result[result$Species == "Empis tessellata", ])
print(result_Empis_tessellata, n = 20)
Empis_tessellata$Empis_tessellata <- c(1,0,5,0,0,0,0,2,0,0,0,3,0,0,1,2,3,0,0,8)

Eristalis_sp <- data.frame(Site = 1:20, Eristalis_sp = rep(0,20))
result_Eristalis<- na.omit(result[result$Species == "Eristalis sp", ])
print(result_Eristalis, n = 20)
Eristalis_sp$Eristalis_sp <- c(2,10,1,2,1,8,0,4,0,1,5,1,0,0,5,1,4,7,1,3)

Eupeodes_corollae <- data.frame(Site = 1:20, Eupeodes_corollae = rep(0,20))
result_Eupeodes_corollae<- na.omit(result[result$Species == "Eupeodes corollae", ])
print(result_Eupeodes_corollae, n = 20)
Eupeodes_corollae$Eupeodes_corollae <- c(3,1,2,4,4,0,0,4,1,2,5,1,1,1,3,9,3,3,2,1)

Helophilus_pendulus <- data.frame(Site = 1:20, Helophilus_pendulus = rep(0,20))
result_Helophilus_pendulus<- na.omit(result[result$Species == "Helophilus pendulus", ])
print(result_Helophilus_pendulus, n = 20)
Helophilus_pendulus$Helophilus_pendulus <- c(1,0,0,0,0,1,0,2,0,0,0,1,2,0,7,3,1,1,0,3)

Lasioglossum_sp <- data.frame(Site = 1:20, Lasioglossum_sp = rep(0,20))
result_Lasioglossum_sp<- na.omit(result[result$Species == "Lasioglossum sp", ])
print(result_Lasioglossum_sp, n = 20)
Lasioglossum_sp$Lasioglossum_sp <- c(3,3,2,4,0,3,1,2,1,1,3,1,1,0,0,0,0,1,6,0)

Maniola_jurtina <- data.frame(Site = 1:20, Maniola_jurtina = rep(0,20))
result_Maniola_jurtina<- na.omit(result[result$Species == "Maniola jurtina", ])
print(result_Maniola_jurtina, n = 20)
Maniola_jurtina$Maniola_jurtina <- c(0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,1,0,3,0,0)

Meligethes_sp <- data.frame(Site = 1:20, Meligethes_sp = rep(0,20))
result_Meligethes_sp<- na.omit(result[result$Species == "Meligethes sp", ])
print(result_Meligethes_sp, n = 20)
Meligethes_sp$Meligethes_sp <- c(0,0,0,2,2,4,5,0,4,0,0,0,6,8,0,3,0,0,0,0)

Merodon_equestris <- data.frame(Site = 1:20, Merodon_equestris = rep(0,20))
result_Merodon_equestris<- na.omit(result[result$Species == "Merodon equestris", ])
print(result_Merodon_equestris, n = 20)
Merodon_equestris$Merodon_equestris <- c(0,2,0,0,0,2,0,0,0,0,1,0,0,0,2,0,0,2,1,0)

Nomada_sp <- data.frame(Site = 1:20, Nomada_sp = rep(0,20))
result_Nomada_sp<- na.omit(result[result$Species == "Nomada sp", ])
print(result_Nomada_sp, n = 20)
Nomada_sp$Nomada_sp <- c(0,1,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0)

Ochlodes_sylvanus <- data.frame(Site = 1:20, Ochlodes_sylvanus = rep(0,20))
result_Ochlodes_sylvanus<- na.omit(result[result$Species == "Ochlodes sylvanus", ])
print(result_Ochlodes_sylvanus, n = 20)
Ochlodes_sylvanus$Ochlodes_sylvanus <- c(0,0,0,1,1,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0)

Oedemera_sp <- data.frame(Site = 1:20, Oedemera_sp = rep(0,20))
result_Oedemera_sp<- na.omit(result[result$Species == "Oedemera sp", ])
print(result_Oedemera_sp, n = 20)
Oedemera_sp$Oedemera_sp <- c(0,0,0,0,6,0,1,0,1,2,0,0,0,0,0,0,0,0,0,0)

Phyllopertha_horticola <- data.frame(Site = 1:20, Phyllopertha_horticola = rep(0,20))
result_Phyllopertha_horticola<- na.omit(result[result$Species == "Phyllopertha horticola", ])
print(result_Phyllopertha_horticola, n = 20)
Phyllopertha_horticola$Phyllopertha_horticola <- c(0,1,0,5,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1)

Sphaerophoria_sp <- data.frame(Site = 1:20, Sphaerophoria_sp = rep(0,20))
result_Sphaerophoria_sp<- na.omit(result[result$Species == "Sphaerophoria sp", ])
print(result_Sphaerophoria_sp, n = 20)
Sphaerophoria_sp$Sphaerophoria_sp <- c(0,2,0,0,2,0,2,0,0,4,0,0,0,0,2,0,0,0,1,0)

Stenurella_melanura <- data.frame(Site = 1:20, Stenurella_melanura = rep(0,20))
result_Stenurella_melanura<- na.omit(result[result$Species == "Stenurella melanura", ])
print(result_Stenurella_melanura, n = 20)
Stenurella_melanura$Stenurella_melanura <- c(5,0,2,2,0,0,6,3,1,0,0,6,1,0,2,2,0,0,2,1)

for(df in list(Andrena_sp, Apis_mellifera, Bombus_pascuorum, Bombus_ruderarius, Bombus_terrestris, 
               Coenonympha_pamphilus, Dasytes_niger, Empis_livida, Empis_tessellata, Eristalis_sp, 
               Eupeodes_corollae, Helophilus_pendulus, Lasioglossum_sp, Maniola_jurtina, Meligethes_sp, 
               Merodon_equestris, Nomada_sp, Ochlodes_sylvanus, Oedemera_sp, Phyllopertha_horticola, 
               Sphaerophoria_sp, Stenurella_melanura)) {
  comb1 <- inner_join(comb1, df, by = "Site")
}

comb1$Arnica_associated <- comb1$Dasytes_niger + comb1$Empis_livida + comb1$Empis_tessellata +
  comb1$Eristalis_sp + comb1$Helophilus_pendulus + comb1$Meligethes_sp + comb1$Merodon_equestris + 
  comb1$Oedemera_sp + comb1$Stenurella_melanura

comb1$Area_associated  <- comb1$Apis_mellifera + comb1$Bombus_ruderarius + comb1$Bombus_terrestris + 
  comb1$Eupeodes_corollae + comb1$Lasioglossum_sp + comb1$Phyllopertha_horticola + 
  comb1$Sphaerophoria_sp + comb1$Andrena_sp + comb1$Bombus_pascuorum + comb1$Coenonympha_pamphilus +
  comb1$Maniola_jurtina + comb1$Ochlodes_sylvanus

comb1$Arnica_associated_narrow <- comb1$Eristalis_sp + comb1$Meligethes_sp +
  comb1$Dasytes_niger + comb1$Empis_livida

comb1$Area_associated_narrow <- comb1$Eupeodes_corollae + comb1$Lasioglossum_sp +
  comb1$Bombus_terrestris + comb1$Bombus_ruderarius + comb1$Bombus_pascuorum +
  comb1$Coenonympha_pamphilus + comb1$Ochlodes_sylvanus

comb1$Arnica_associated_noMel <- comb1$Dasytes_niger + comb1$Empis_livida + comb1$Empis_tessellata +
  comb1$Eristalis_sp + comb1$Helophilus_pendulus + comb1$Merodon_equestris + 
  comb1$Oedemera_sp + comb1$Stenurella_melanura

#add number of Arnica by multiplying with overall number of pollen
comb_all$Nr_Arnica <- comb_all$nPoll*comb_all$P_ASTE.Arnica.montana
comb_all$Nr_Not.Arnica <- comb_all$nPoll*comb_all$P_Not.Arnica

#list of species with 5 or more samples
work_species <- as.list(imp_species$Species)
work_species

#second combined dataframe with only those important species
comb_all2 <- comb_all[comb_all$Species %in% work_species, ]
which(!complete.cases(comb_all2)) #find indices of rows with NA's
#comb_all2 <- comb_all2[-373,] #remove Andrena sample with no pollen identified
comb_all2$Nr_Arnica <- round(comb_all2$Nr_Arnica)
comb_all2$Nr_Not.Arnica <- round(comb_all2$Nr_Not.Arnica)
comb_all2$nPoll <- round(comb_all2$nPoll)
comb_all2$nGarb <- round(comb_all2$nGarb)

#remove unnecessary objects
rm(adj)
rm(Andrena_sp)
rm(Apis_mellifera)
rm(av_seeds)
rm(Bombus_pascuorum)
rm(Bombus_ruderarius)
rm(Bombus_terrestris)
rm(Coenonympha_pamphilus)
rm(Dasytes_niger)
rm(df)
rm(Empis_livida)
rm(Empis_tessellata)
rm(Eristalis_sp)
rm(Eupeodes_corollae)
rm(Helophilus_pendulus)
rm(Lasioglossum_sp)
rm(Maniola_jurtina)
rm(Meligethes_sp)
rm(Merodon_equestris)
rm(Nomada_sp)
rm(Ochlodes_sylvanus)
rm(Oedemera_sp)
rm(ori)
rm(Phyllopertha_horticola)
rm(polli)
rm(result)
rm(result_Andrena)
rm(result_Apis_mellifera)
rm(result_Bombus_pascuorum)
rm(result_Bombus_ruderarius)
rm(result_Bombus_terrestris)
rm(result_Coenonympha_pamphilus)
rm(result_Dasytes_niger)
rm(result_Empis_livida)
rm(result_Empis_tessellata)
rm(result_Eristalis)
rm(result_Eupeodes_corollae)
rm(result_Helophilus_pendulus)
rm(result_Lasioglossum_sp)
rm(result_Maniola_jurtina)
rm(result_Meligethes_sp)
rm(result_Merodon_equestris)
rm(result_Nomada_sp)
rm(result_Ochlodes_sylvanus)
rm(result_Oedemera_sp)
rm(result_Phyllopertha_horticola)
rm(result_Sphaerophoria_sp)
rm(result_Stenurella_melanura)
rm(seeds)
rm(species_counts)
rm(Sphaerophoria_sp)
rm(Stenurella_melanura)
rm(work_species)
gc()

#subset into pollinators associated with Arnica and pollinators rather associated with the area
Arnica_polli <- subset(comb_all2, Species %in% c("Eristalis sp", "Meligethes sp", 
                                                 "Stenurella melanura", "Dasytes niger", 
                                                 "Empis livida", "Empis tessellata", 
                                                 "Helophilus pendulus", "Oedemera sp",
                                                 "Merodon equestris"))

Area_polli <- subset(comb_all2, Species %in% c("Bombus terrestris", "Bombus ruderarius", 
                                               "Lasioglossum sp", "Apis mellifera",
                                               "Eupeodes corollae", "Sphaerophoria sp",
                                               "Andrena sp", "Phyllopertha horticola",
                                               "Bombus pascuorum", "Coenonympha pamphilus",
                                               "Maniola jurtina", "Ochlodes sylvanus"))

#add association to comb_all2
comb_all2$Association <- rep(0,length(comb_all2$Site))

for (i in 1:length(comb_all2$Site)){
  comb_all2[i, "Association"] <- ifelse(comb_all2[i,"Species"] %in% c("Eristalis sp", "Meligethes sp", 
                                                                          "Stenurella melanura", "Dasytes niger", 
                                                                          "Empis livida", "Empis tessellata", 
                                                                          "Helophilus pendulus", "Oedemera sp",
                                                                          "Merodon equestris"),
                                            "Arnica", "Area")
}


#load meteorological data
meteo <- read.csv("meteo_extended.csv", h = T) #load meteorological observations and  # pollinators per site
meteo$Site <- meteo$Site_id

poll_nr <- inner_join(pop, meteo, by="Site")
poll_nr <- poll_nr[,-11]

#sampling time as qualitative:
poll_nr$time_quali <- rep(NA, length(poll_nr$Visit_id))
poll_nr[poll_nr$End_time_R < 11, ]$time_quali <- "morning"
poll_nr[poll_nr$End_time_R >= 11 & poll_nr$End_time_R < 13, ]$time_quali <- "early_midday"
poll_nr[poll_nr$End_time_R >= 13 & poll_nr$End_time_R < 15, ]$time_quali <- "late_midday"
poll_nr[poll_nr$End_time_R >= 15 & poll_nr$End_time_R < 17, ]$time_quali <- "afternoon"
poll_nr[poll_nr$End_time_R >= 17, ]$time_quali <- "evening"

rm(meteo)
rm(pop)
rm(comb2)
gc()
