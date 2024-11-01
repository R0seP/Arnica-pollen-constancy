#This script aims to investigate the pollinator-plant network of the pollinators
#caught at Arnica sites using bipartite graphs.

library(tidyverse)
library(vegan)
library(statnet.common)
library(network)
library(sna)
library(bipartite)

#prepare data----
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)

network_data <- comb_all2 #create new data frame for this analysis

for (r in 1:length(network_data$Site)){
  for (c in 1:34){
    network_data[r,14+c] <- comb_all2[r,14+c]*comb_all2[r, "nPoll"]
  }
} #calculate number of pollen per pollen group

network_data <- network_data[,-c(1:3,5:14,50:53)]  #delete not needed columns
network_data <- na.omit(network_data) #remove NAs

network_sums <- network_data %>%
  group_by(Species) %>%
  summarize(across(everything(), sum)) #sum over the number of pollen of each plant per pollinator species

network_sums[,c(2:36)] <- round(network_sums[,c(2:36)]) #round to an integer number of pollen
names(network_sums) <- c("Pollinator", "Viburnum", "Apiaceae", "Arnica_montana", 
                         "Cirsium", "Crepis", "Matricaria", "Senecio", "Echium_vulgare",
                         "Brassicaceae", "Lonicera_periclymenum", "Succisa_pratensis",
                         "Caryophyllaceae", "Ericaceae", "Astralagus", "Lathyrus",
                         "Trifolium_pratense", "Trifolium_repens", "Vicia", "Geranium",
                         "Prunella", "Tilia", "Epilobium", "Rhinanthus", "Digitalis",
                         "Veronica", "Poaceae", "Rumex", "Ranunculus", "Filipendula",
                         "Geum", "Potentilla", "Rosa", "Sorbus", "Spirea", "Galium")

network <- as.matrix(network_sums[,2:36]) #turn network data into matrix for further use
rownames(network) <- network_sums$Pollinator

#visualize networks----
visweb(network)
plotweb(network, text.rot=90, col.low = "steelblue", col.high = "forestgreen")
