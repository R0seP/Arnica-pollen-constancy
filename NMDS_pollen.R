#This script aims to look at the differences in pollen community composition
#between pollinators of different species

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")

library(vegan)
library(tidyverse)

#get data----
source("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R/data_preparation.R", echo = TRUE)
adj <- read.csv("adj_new.csv", h = T)
polli <- read.csv("pollinators.csv", h = T) 

nmds_pollen <- adj #create new data frame for this analysis

for (r in 1:length(nmds_pollen$Site)){
  for (c in 1:34){
    nmds_pollen[r,7+c] <- adj[r,7+c]*adj[r, "nPoll"]
  }
} #calculate number of pollen per pollen group

nmds_pollen <- inner_join(nmds_pollen, polli, by = c("Site", "Pollinator")) #add pollinator species info
nmds_pollen <- nmds_pollen[,-c(2:7,46,48)]  #delete not needed columns
nmds_pollen <- nmds_pollen[,c(1,37:40,2:36)] #resort 
nmds_pollen <- na.omit(nmds_pollen)

work_species <- as.list(imp_species$Species)
nmds_pollen <- nmds_pollen[nmds_pollen$Species %in% work_species, ]
#exclude all but the 22 "important species" with more than 5 obseravtions

#NMDS----
NMDS <- metaMDS(as.matrix(nmds_pollen[,6:40]), distance = "bray", k = 3, autotransform = TRUE, trymax=100)
NMDS$stress
#reliable NMDS if stress < 0.2, stress currently approximately 0.2 -> still ok,
#especially because so much data?

# make new dataframe with by extracting NMDS scores
nmds.scores <- as.data.frame(scores(NMDS)$sites)
nmds.scores <- nmds.scores %>%
  mutate(Site = as.factor(nmds_pollen$Site), Group = as.factor(nmds_pollen$Group),
         Species = as.factor(nmds_pollen$Species))

#visualization----
# define hidden vegan function that finds coordinates for drawing a covariance ellipse
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


# create empty dataframe to combine NMDS data with ellipse data
ellipse_df <- data.frame()

# adding data for ellipse, in this case using distance as a grouping factor
for(g in levels(nmds.scores$Species)){
  ellipse_df <- rbind(ellipse_df, cbind(as.data.frame(with(nmds.scores[nmds.scores$Species==g,],
                                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                  wt=rep(1/length(NMDS1),length(NMDS1)))$cov,
                                                                           center=c(mean(NMDS1),mean(NMDS2)))))
                                        ,Distance=g))
}


# create ggplot
red_palette <- c("#FF0000", "#E60000", "#CC0000", "#B30000", "#990000", "#800000", "#660000") #bees
green_palette <- c("#66C2A5", "#1B9E77", "#00441B") #butterflies
blue_palette <- c("#1F77B4", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7") #beetles
yellow_palette <- c("#FFFF00", "#FFEA00", "#FFD700", "#FFC300", "#FFB000", "#FF9C00", "#FF8800") #flies


color_palette <- c("#FF0000", "#E60000", "#CC0000", "#B30000", "#990000", "#66C2A5", 
                   "#1F77B4", "#FFFF00", "#FFEA00", "#FFD700", "#FFC300", "#FFB000", 
                   "#800000", "#1B9E77", "#6BAED6", "#FF9C00", "#660000", "#00441B", 
                   "#9ECAE1", "#C6DBEF", "#FF8800", "#DEEBF7")

(NMDS_plot <- ggplot(data = nmds.scores, aes(NMDS1, NMDS2)) +
    geom_point(aes(color = Species, shape = Species)) + # adding different colors and shapes for different groups
    geom_path(data=ellipse_df, aes(x=NMDS1, y=NMDS2, colour=Distance), linewidth=1) + # adding covariance ellipses according to group
    scale_shape_manual(values = 1:22) + # manually specifying shapes for 22 species
    scale_color_manual(values = color_palette) + # manually specifying colors for 22 species
    guides(color = guide_legend(override.aes = list(linetype=c(NA,NA)))) + # removes lines from legend
    theme_bw() + # adding theme
    theme(panel.grid = element_blank(), # remove background grid
          legend.text = element_text(size = 12), # increase legend text size
          legend.title = element_text(size = 14), # increase legend title size
          axis.text = element_text(size = 12), # increase axis text size
          axis.title = element_text(size = 14)))  # increase axis title size) 


#permanova----
# using a PERMANOVA (PERmutational Multivariate ANalysis Of VAriance) 
# to test the differences in community composition
permanova_species <- adonis2(as.matrix(nmds_pollen [,6:40]) ~ Species + Group + Site
                             , nmds_pollen, permutations = 999, method = "bray") 
# permutations is the number of possible arrangements the data could take
# using permutations = 999 is standard practice in science and across the literature
# can increase to get a more thorough analysis
# use method = bray as it is what we previously used to calculate pairwise distances

permanova_species


# check model assumptions
# check for multivariate homogeneity of group variances
# generate distance matrix from pollinator community matrix
dist <- vegdist(as.matrix(nmds_pollen [,6:40]), method = "bray")

# use betadisper test to check for multivariate homogeneity of group variances
dispersion <- betadisper(dist, group=nmds_pollen$Species)
permutest(dispersion)
# test gives a significant result, meaning that variances are not homogeneous,
# model might not be reliable

# Extract distances to centroids
dispersion_distances <- dispersion$distances

# Add dispersion distances to data
nmds_pollen$Dispersion <- dispersion_distances

#add dispersion to model to account for non-homogeneous variances
permanova_species <- adonis2(as.matrix(nmds_pollen [,6:40]) ~ Species + Group + Site + Dispersion
                             , nmds_pollen, permutations = 999, method = "bray") 
permanova_species

#simper----
# simper (SIMilarity PERcentage) analysis
# compare community differences, report what species are driving those differences.
nmds_matrix <- as.matrix(nmds_pollen[, 6:40])
simper <- with(nmds_pollen[,c(1:5)], simper(nmds_matrix), Species, permutations = 100)

# see most influential species contributing to community differences
simper
#Species represented in this output are responsible for 70% of the observed variation. 
#Each time it adds another taxonomic group it shows its cumulative contribution of 
#rather than its individual one.

summary(simper)
