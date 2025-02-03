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
comb_all2 <- read.csv("comb_all2.csv", h = T)

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
network <- t(network) #transmute matrix to have higher trophic level as column names

###visualize network----
visweb(network, labsize = 1.5)
plotweb(network, text.rot=90, col.low = "forestgreen", col.high = "steelblue",
        labsize = 1.5)


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
  site_matrix <- t(site_matrix) #transpose to make bipartite recognize higher trophic level
  
  return(site_matrix)
}

#apply the function to all sites
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
plotweb(web1, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#2
visweb(web2)
plotweb(web2, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#3
visweb(web3)
plotweb(web3, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#4
visweb(web4)
plotweb(web4, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#5
visweb(web5)
plotweb(web5, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#6
visweb(web6)
plotweb(web6, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#7
visweb(web7)
plotweb(web7, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#8
visweb(web8)
plotweb(web8, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#9
visweb(web9)
plotweb(web9, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#10
visweb(web10)
plotweb(web10, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#11
visweb(web11)
plotweb(web11, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#12
visweb(web12)
plotweb(web12, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#13
visweb(web13)
plotweb(web13, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#14
visweb(web14)
plotweb(web14, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#15
visweb(web15)
plotweb(web15, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#16
visweb(web16)
plotweb(web16, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#17
visweb(web17)
plotweb(web17, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#18
visweb(web18)
plotweb(web18, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#19
visweb(web19)
plotweb(web19, text.rot=90, col.low = "forestgreen", col.high = "steelblue")

#20
visweb(web20)
plotweb(web20, text.rot=90, col.low = "forestgreen", col.high = "steelblue")


###species-level metrics----
species_metrics <- specieslevel(network)
species_metrics

###calculate Blüthgen's d
#(specialization index) for Arnica at different sites

d_webs <- lapply(webs, specieslevel, index = 'd')

#vector to store d values
d_values <- numeric(length(d_webs))

#extract d-value for Arnica across all sites
for (i in seq_along(d_webs)) {
  site <- d_webs[[i]]
  d_values[i] <- site$`lower level`[plant_species, "d"]
}

d_values

#null models for all webs
set.seed(2800)  # For reproducibility
null_models1 <- replicate(1000, r2dtable(1, rowSums(web1), colSums(web1))[[1]])
null_models2 <- replicate(1000, r2dtable(1, rowSums(web2), colSums(web2))[[1]])
null_models3 <- replicate(1000, r2dtable(1, rowSums(web3), colSums(web3))[[1]])
null_models4 <- replicate(1000, r2dtable(1, rowSums(web4), colSums(web4))[[1]])
null_models5 <- replicate(1000, r2dtable(1, rowSums(web5), colSums(web5))[[1]])
null_models6 <- replicate(1000, r2dtable(1, rowSums(web6), colSums(web6))[[1]])
null_models7 <- replicate(1000, r2dtable(1, rowSums(web7), colSums(web7))[[1]])
null_models8 <- replicate(1000, r2dtable(1, rowSums(web8), colSums(web8))[[1]])
null_models9 <- replicate(1000, r2dtable(1, rowSums(web9), colSums(web9))[[1]])
null_models10 <- replicate(1000, r2dtable(1, rowSums(web10), colSums(web10))[[1]])
null_models11 <- replicate(1000, r2dtable(1, rowSums(web11), colSums(web11))[[1]])
null_models12 <- replicate(1000, r2dtable(1, rowSums(web12), colSums(web12))[[1]])
null_models13 <- replicate(1000, r2dtable(1, rowSums(web13), colSums(web13))[[1]])
null_models14 <- replicate(1000, r2dtable(1, rowSums(web14), colSums(web14))[[1]])
null_models15 <- replicate(1000, r2dtable(1, rowSums(web15), colSums(web15))[[1]])
null_models16 <- replicate(1000, r2dtable(1, rowSums(web16), colSums(web16))[[1]])
null_models17 <- replicate(1000, r2dtable(1, rowSums(web17), colSums(web17))[[1]])
null_models18 <- replicate(1000, r2dtable(1, rowSums(web18), colSums(web18))[[1]])
null_models19 <- replicate(1000, r2dtable(1, rowSums(web19), colSums(web19))[[1]])
null_models20 <- replicate(1000, r2dtable(1, rowSums(web20), colSums(web20))[[1]])

#calculate Blüthgen's d for each null model
null_d_1 <- apply(null_models1, 3, function(x) specieslevel(x, index = "d"))
null_d_2 <- apply(null_models2, 3, function(x) specieslevel(x, index = "d"))
null_d_3 <- apply(null_models3, 3, function(x) specieslevel(x, index = "d"))
null_d_4 <- apply(null_models4, 3, function(x) specieslevel(x, index = "d"))
null_d_5 <- apply(null_models5, 3, function(x) specieslevel(x, index = "d"))
null_d_6 <- apply(null_models6, 3, function(x) specieslevel(x, index = "d"))
null_d_7 <- apply(null_models7, 3, function(x) specieslevel(x, index = "d"))
null_d_8 <- apply(null_models8, 3, function(x) specieslevel(x, index = "d"))
null_d_9 <- apply(null_models9, 3, function(x) specieslevel(x, index = "d"))
null_d_10 <- apply(null_models10, 3, function(x) specieslevel(x, index = "d"))
null_d_11 <- apply(null_models11, 3, function(x) specieslevel(x, index = "d"))
null_d_12 <- apply(null_models12, 3, function(x) specieslevel(x, index = "d"))
null_d_13 <- apply(null_models13, 3, function(x) specieslevel(x, index = "d"))
null_d_14 <- apply(null_models14, 3, function(x) specieslevel(x, index = "d"))
null_d_15 <- apply(null_models15, 3, function(x) specieslevel(x, index = "d"))
null_d_16 <- apply(null_models16, 3, function(x) specieslevel(x, index = "d"))
null_d_17 <- apply(null_models17, 3, function(x) specieslevel(x, index = "d"))
null_d_18 <- apply(null_models18, 3, function(x) specieslevel(x, index = "d"))
null_d_19 <- apply(null_models19, 3, function(x) specieslevel(x, index = "d"))
null_d_20 <- apply(null_models20, 3, function(x) specieslevel(x, index = "d"))

#extract d-values only for Arnica
Arnica_null_d1 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d1[i] <- null_d_1[[i]][[2]]$d[3]
}

Arnica_null_d2 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d2[i] <- null_d_2[[i]][[2]]$d[3]
}

Arnica_null_d3 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d3[i] <- null_d_3[[i]][[2]]$d[3]
}

Arnica_null_d4 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d4[i] <- null_d_4[[i]][[2]]$d[3]
}

Arnica_null_d5 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d5[i] <- null_d_5[[i]][[2]]$d[3]
}

Arnica_null_d6 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d6[i] <- null_d_6[[i]][[2]]$d[3]
}

Arnica_null_d7 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d7[i] <- null_d_7[[i]][[2]]$d[3]
}

Arnica_null_d8 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d8[i] <- null_d_8[[i]][[2]]$d[3]
}

Arnica_null_d9 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d9[i] <- null_d_9[[i]][[2]]$d[3]
}

Arnica_null_d10 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d10[i] <- null_d_10[[i]][[2]]$d[3]
}

Arnica_null_d11 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d11[i] <- null_d_11[[i]][[2]]$d[3]
}

Arnica_null_d12 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d12[i] <- null_d_12[[i]][[2]]$d[3]
}

Arnica_null_d13 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d13[i] <- null_d_13[[i]][[2]]$d[3]
}

Arnica_null_d14 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d14[i] <- null_d_14[[i]][[2]]$d[3]
}

Arnica_null_d15 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d15[i] <- null_d_15[[i]][[2]]$d[3]
}

Arnica_null_d16 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d16[i] <- null_d_16[[i]][[2]]$d[3]
}

Arnica_null_d17 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d17[i] <- null_d_17[[i]][[2]]$d[3]
}

Arnica_null_d18 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d18[i] <- null_d_18[[i]][[2]]$d[3]
}

Arnica_null_d19 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d19[i] <- null_d_19[[i]][[2]]$d[3]
}

Arnica_null_d20 <- rep(0, 1000)
for(i in 1:1000){
  Arnica_null_d20[i] <- null_d_20[[i]][[2]]$d[3]
}

#calculate delta_d, delta transformation by standardizing d against null models
delta_d1 <- mean(Arnica_null_d1) - d_values[1]
delta_d2 <- mean(Arnica_null_d2) - d_values[2]
delta_d3 <- mean(Arnica_null_d3) - d_values[3]
delta_d4 <- mean(Arnica_null_d4) - d_values[4]
delta_d5 <- mean(Arnica_null_d5) - d_values[5]
delta_d6 <- mean(Arnica_null_d6) - d_values[6]
delta_d7 <- mean(Arnica_null_d7) - d_values[7]
delta_d8 <- mean(Arnica_null_d8) - d_values[8]
delta_d9 <- mean(Arnica_null_d9) - d_values[9]
delta_d10 <- mean(Arnica_null_d10) - d_values[10]
delta_d11 <- mean(Arnica_null_d11) - d_values[11]
delta_d12 <- mean(Arnica_null_d12) - d_values[12]
delta_d13 <- mean(Arnica_null_d13) - d_values[13]
delta_d14 <- mean(Arnica_null_d14) - d_values[14]
delta_d15 <- mean(Arnica_null_d15) - d_values[15]
delta_d16 <- mean(Arnica_null_d16) - d_values[16]
delta_d17 <- mean(Arnica_null_d17) - d_values[17]
delta_d18 <- mean(Arnica_null_d18) - d_values[18]
delta_d19 <- mean(Arnica_null_d19) - d_values[19]
delta_d20 <- mean(Arnica_null_d20) - d_values[20]

delta_d <- list(delta_d1, delta_d2, delta_d3, delta_d4, delta_d5, delta_d6, delta_d7,
                delta_d8, delta_d9, delta_d10, delta_d11, delta_d12, delta_d13, delta_d14,
                delta_d15, delta_d16, delta_d17, delta_d18, delta_d19, delta_d20)

###combine species-level metrics in data frame.
comb1 <- read.csv("comb1.csv", h = T)
species_metrics <- comb1[,1:2]
species_metrics$Arnica_d <- as.numeric(d_values)
species_metrics$Arnica_delta_d <- abs(as.numeric(delta_d))
#by taking the absolute value of delta d webs that had high d-values now again 
#have high d-values because null d's were very small

#save species-level metrics for further analysis
write.csv(species_metrics, "species_metrics.csv", row.names = F)
