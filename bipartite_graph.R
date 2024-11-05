#This script aims to investigate the pollinator-plant network of the pollinators
#caught at Arnica sites using bipartite graphs.

library(tidyverse)
library(vegan)
library(statnet.common)
library(network)
library(sna)
library(bipartite)

#entire network----
###prepare data----
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

###visualize network----
visweb(network)
plotweb(network, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

###calculate network metrics----
#generate null models
set.seed(2800)  # For reproducibility
null_models <- replicate(1000, r2dtable(1, rowSums(network), colSums(network))[[1]])

#calculate nestedness
nestedness <- networklevel(network, index = "NODF")
#nestedness 78 out of 100, relatively high?

#calculate nestedness for each null model
null_nestedness <- apply(null_models, 3, function(x) networklevel(x, index = "NODF"))

#perform t-test for nestedness
t_test_nest <- t.test(null_nestedness, mu = nestedness)
t_test_nest
#significantly nested network (significantly higher NODF than expected by chance)

#calculate links per species
links_per_species <- networklevel(network, index = "links per species")

#calculate links per species for each null model
null_links <- apply(null_models, 3, function(x) networklevel(x, index = "links per species"))

#perform t-test for links per species
t_test_links <- t.test(null_links, mu = links_per_species)
t_test_links
#significantly fewer links per species in network than expected by chance


#calculate H2 (specialization)
#calculate H2
H2 <- networklevel(network, index = "H2")

#calculate H2 for each null model
null_H2 <- apply(null_models, 3, function(x) networklevel(x, index = "H2"))

#perform t-test for H2
t_test_H2 <- t.test(null_H2, mu = H2)
t_test_H2
#significantly higher specialization in observed network than expected by chance


#site networks----
###prepare data----
network_data <- comb_all2 #create new data frame for this analysis

for (r in 1:length(network_data$Site)){
  for (c in 1:34){
    network_data[r,14+c] <- comb_all2[r,14+c]*comb_all2[r, "nPoll"]
  }
} #calculate number of pollen per pollen group

network_data <- network_data[,-c(2:3,5:14,50:53)]  #delete not needed columns
network_data <- na.omit(network_data) #remove NAs

#network data per site
#create function to prepare data for each site
process_site_data <- function(site_number, network_data) {
  site_data <- subset(network_data, Site == site_number)
  site_data <- site_data[,-1]
  
  site_sums <- site_data %>%
    group_by(Species) %>%
    summarize(across(everything(), sum)) # Sum over the number of pollen of each plant per pollinator species
  
  site_sums[,c(2:36)] <- round(site_sums[,c(2:36)]) # Round to an integer number of pollen
  names(site_sums) <- c("Pollinator", "Viburnum", "Apiaceae", "Arnica_montana", 
                        "Cirsium", "Crepis", "Matricaria", "Senecio", "Echium_vulgare",
                        "Brassicaceae", "Lonicera_periclymenum", "Succisa_pratensis",
                        "Caryophyllaceae", "Ericaceae", "Astralagus", "Lathyrus",
                        "Trifolium_pratense", "Trifolium_repens", "Vicia", "Geranium",
                        "Prunella", "Tilia", "Epilobium", "Rhinanthus", "Digitalis",
                        "Veronica", "Poaceae", "Rumex", "Ranunculus", "Filipendula",
                        "Geum", "Potentilla", "Rosa", "Sorbus", "Spirea", "Galium")
  
  site_matrix <- as.matrix(site_sums[,2:36]) # Turn network data into matrix for further use
  rownames(site_matrix) <- site_sums$Pollinator
  
  return(site_matrix)
}

#apply the function to all sites (assuming you have 20 sites)
site_numbers <- 1:20
site_matrices <- lapply(site_numbers, process_site_data, network_data = network_data)

#create matrices for each site
web1 <- site_matrices[[1]]
web2 <- site_matrices[[2]]
web3 <- site_matrices[[3]]
web4 <- site_matrices[[4]]
web5 <- site_matrices[[5]]
web6 <- site_matrices[[6]]
web7 <- site_matrices[[7]]
web8 <- site_matrices[[8]]
web9 <- site_matrices[[9]]
web10 <- site_matrices[[10]]
web11 <- site_matrices[[11]]
web12 <- site_matrices[[12]]
web13 <- site_matrices[[13]]
web14 <- site_matrices[[14]]
web15 <- site_matrices[[15]]
web16 <- site_matrices[[16]]
web17 <- site_matrices[[17]]
web18 <- site_matrices[[18]]
web19 <- site_matrices[[19]]
web20 <- site_matrices[[20]]

webs <- list(web1, web2, web3, web4, web5, web6, web7, web8, web9, web10, 
             web11, web12, web13, web14, web15, web16, web17, web18, web19, web20)


###visualize networks----
#1
visweb(web1)
plotweb(web1, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#2
visweb(web2)
plotweb(web2, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#3
visweb(web3)
plotweb(web3, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#4
visweb(web4)
plotweb(web4, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#5
visweb(web5)
plotweb(web5, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#6
visweb(web6)
plotweb(web6, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#7
visweb(web7)
plotweb(web7, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#8
visweb(web8)
plotweb(web8, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#9
visweb(web9)
plotweb(web9, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#10
visweb(web10)
plotweb(web10, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#11
visweb(web11)
plotweb(web11, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#12
visweb(web12)
plotweb(web12, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#13
visweb(web13)
plotweb(web13, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#14
visweb(web14)
plotweb(web14, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#15
visweb(web15)
plotweb(web15, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#16
visweb(web16)
plotweb(web16, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#17
visweb(web17)
plotweb(web17, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#18
visweb(web18)
plotweb(web18, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#19
visweb(web19)
plotweb(web19, text.rot=90, col.low = "steelblue", col.high = "forestgreen")

#20
visweb(web20)
plotweb(web20, text.rot=90, col.low = "steelblue", col.high = "forestgreen")


###calculate network metrics----
nested_webs <- lapply(webs, networklevel, index = 'NODF')
links_webs <- lapply(webs, networklevel, index = 'links per species')
H2_webs <- lapply(webs, networklevel, index = 'H2')

#combine network metrics in data frame.
network_metrics <- comb1[,1:2]
network_metrics$NODF <- as.numeric(nested_webs)
network_metrics$links_per_species <- as.numeric(links_webs)
network_metrics$H2 <- as.numeric(H2_webs)

head(network_metrics)
str(network_metrics)

#save network metrics for further analysis
write.csv(network_metrics, "network_metrics.csv", row.names = F)
