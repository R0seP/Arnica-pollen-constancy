#This script takes the output of the pollen identification with adjusted classifications
#and modifies it so that the garbage category "Clumps" is no longer counted in 
#the number and percentage of pollen grains in a sample. 
#original data with adjusted classifications
adj <- read.csv("pollen_adj_class_without_Senecio.csv", h = T)

#create new data frame 
adj_new <- adj

#correct the number of pollen by subtracting the number of clumps
adj_new$nPoll <- adj$nPoll - (adj$nPoll * adj$P_zzClumps)
#correct number of garbage objects by adding the number of clumps
adj_new$nGarb <- adj$nGarb + (adj$nPoll * adj$P_zzClumps)

#add column to old data frame that sums up proportions of all pollen (but not the clumps)
adj$P_pollen <- rowSums(adj[, 8:41], na.rm = TRUE)

#adjust the percentages of pollen grains so that they add up to 100% again
for (r in 1:length(adj_new$Site)){
  for (c in 1:35){
    adj_new[r,7+c] <- adj[r,7+c]/adj[r,"P_pollen"]
  }
}

#check if loop worked and all row sums in new data frame add up to 1 without clumps
rowSums(adj_new[, 8:41], na.rm = TRUE)
  #add to 1 but some 0's there because of no pollen data
  #okay because I want to keep original data, but note which ones:
which(rowSums(adj_new[, 8:42], na.rm = TRUE) == 0)
#row 216, row 296, row 352 and row 468 contain NA's because of no pollen recognized in the sample

#delete clumps percentages from new data frame
adj_new <- adj_new[,-42]


#save new data frame
write.csv(adj_new, "adj_new_without_Senecio.csv", row.names = F)
