#this script summarizes the number of individuals caught and the percentage of
#Arnica on the 22 pollinator species with 5 or more occurences

library(dplyr)
setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")

source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)

#count the number of individuals per species
##total 
species_counts_total <- comb_all2%>%
  group_by(Species)%>%
  summarise(Count = n())

summary <- as.data.frame(species_counts_total)

names(summary) <- c("Species", "Nr.Pollinators_total")

##on Arnica
species_counts_Arnica <- subset(comb_all2, Group == "flower")%>%
  group_by(Species)%>%
  summarise(Count = n())

species_counts_Arnica <- as.data.frame(species_counts_Arnica)

names(species_counts_Arnica) <- c("Species", "Nr.Pollinators_on_Arnica")

summary <- left_join(summary, species_counts_Arnica, by = "Species")

##in area
species_counts_area <- subset(comb_all2, Group == "area")%>%
  group_by(Species)%>%
  summarise(Count = n())

species_counts_area <- as.data.frame(species_counts_area)

names(species_counts_area) <- c("Species", "Nr.Pollinators_in_area")

summary <- left_join(summary, species_counts_area, by = "Species")



#calculate means for insects caught on Arnica and in area separately
##total
species_means <- comb_all2 %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE))

species_means <- as.data.frame(species_means)
names(species_means) <- c("Species", "Perc.Arnica_total")

summary <- left_join(summary, species_means, by = "Species")

##Arnica-caught
species_means_Arnica <- subset(comb_all2, Group == "flower") %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE))

species_means_Arnica <- as.data.frame(species_means_Arnica)
names(species_means_Arnica) <- c("Species", "Perc.Arnica_on_Arnica")

summary <- left_join(summary, species_means_Arnica, by = "Species")

##Area-caught
species_means_area <- subset(comb_all2, Group == "area") %>%
  group_by(Species) %>%
  summarize(mean_value = mean(P_ASTE.Arnica.montana, na.rm = TRUE))

species_means_area <- as.data.frame(species_means_area)
names(species_means_area) <- c("Species", "Perc.Arnica_in_area")

summary <- left_join(summary, species_means_area, by = "Species")


#add row with totals
sum_row <- colSums(summary[,2:4], na.rm = TRUE)

overall_mean <- mean(summary$Perc.Arnica_total)
Arnica_mean <- mean(summary$Perc.Arnica_on_Arnica, na.rm = TRUE)
Area_mean <- mean(summary$Perc.Arnica_in_area, na.rm = TRUE)

sum_row <- c("total", sum_row, overall_mean, Arnica_mean, Area_mean)

summary <- rbind(summary, sum_row)
summary$Nr.Pollinators_total <- as.numeric(summary$Nr.Pollinators_total)
summary$Nr.Pollinators_on_Arnica <- as.numeric(summary$Nr.Pollinators_on_Arnica)
summary$Nr.Pollinators_in_area <- as.numeric(summary$Nr.Pollinators_in_area)
summary$Perc.Arnica_total <- as.numeric(summary$Perc.Arnica_total)
summary$Perc.Arnica_on_Arnica <- as.numeric(summary$Perc.Arnica_on_Arnica)
summary$Perc.Arnica_in_area <- as.numeric(summary$Perc.Arnica_in_area)

#multiply percentage columns with 100 to get percent
summary$Perc.Arnica_total <- summary$Perc.Arnica_total * 100
summary$Perc.Arnica_on_Arnica <- summary$Perc.Arnica_on_Arnica * 100
summary$Perc.Arnica_in_area <- summary$Perc.Arnica_in_area * 100

#save summary table as .csv
write.csv(summary, "pollinators_summarized.csv", row.names = F)
