##This script aims to analyze the pollinator community composition both in the
##pollinators caught on Arnica and in the pollinators caught in the area

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")

library(vegan)
library(tidyverse)


#visualization----
polli <- read.csv("pollinators.csv", h = T)

par(mfrow = c(1,1))
mosaicplot(Species ~ Group, data = polli,
           las = 2,
           main = "",
           ylab = "Group",
           xlab = "Pollinator species",
           col = c("darkgreen", "orange"),
           cex = 0.9)

mosaicplot(Family ~ Group, data = polli,
           las = 2,
           main = "",
           ylab = "Group",
           xlab = "Pollinator family",
           col = c("darkgreen", "orange"),
           cex = 1.2)

mosaicplot(Overview ~ Group, data = polli,
           las = 2,
           main = "",
           ylab = "Group",
           xlab = "Pollinator class",
           col = c("darkgreen", "orange"),
           cex = 1.2)

#only species with more than 5 occurrences
mosaicplot(Species ~ Group, data = comb_all2,
           las = 2,
           main = "",
           ylab = "Group",
           xlab = "Pollinator species",
           col = c("darkgreen","orange"),
           cex = 1.2)


#NMDS----
#Non-metric multidimensional scaling (NMDS) to compare species composition on 
#Arnica and in area using Bray-Curtis similarity metric

#what influences the pollinator composition?

nmds_data <- read.csv("species_numbers_persite.csv")
community <- nmds_data[,3:85]
community_matrix <- as.matrix(community)
NMDS <- metaMDS(community_matrix, distance = "bray", k = 3, autotransform = TRUE, trymax=100)
NMDS$stress
#stress alright (<0.2) when k = 3, stress a bit high (~ 0.24) when k = 2

#plot(NMDS)
#ordiplot(NMDS, type = "n") # create blank ordination plot
#orditorp(NMDS, display = "species", col="red", air = 0.1) # add species names in red
#orditorp(NMDS, display = "sites", cex = 1.25, air = 0.1) # add site numbers in black

#ordiplot(NMDS) # plot shows communities (circles) and species (crosses)
#ordiellipse(NMDS, nmds_data$Group, label = FALSE,
#            col=c("darkgreen", "orange"),
#            draw = "polygon", alpha=120) # adding ellipses to the plot, grouping by group (nmds_data$Group)
#legend("topright", title="Group",
#       c("area","flower"), fill=c("darkgreen", "orange"), horiz=FALSE, cex=.9) # adding a legend


#difference groups----
# make new dataframe with by extracting NMDS scores
nmds.scores <- as.data.frame(scores(NMDS)$sites)
nmds.scores <- nmds.scores %>%
  mutate(Site = as.factor(nmds_data$Site), Group = as.factor(nmds_data$Group))

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
for(g in levels(nmds.scores$Group)){
  ellipse_df <- rbind(ellipse_df, cbind(as.data.frame(with(nmds.scores[nmds.scores$Group==g,],
                                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                  wt=rep(1/length(NMDS1),length(NMDS1)))$cov,
                                                                           center=c(mean(NMDS1),mean(NMDS2)))))
                                        ,Distance=g))
}


# create ggplot
(NMDS_plot <- ggplot(data = nmds.scores, aes(NMDS1, NMDS2)) +
    geom_point(aes(color = Group, shape = Group)) + # adding different colors and shapes for different groups
    geom_path(data=ellipse_df, aes(x=NMDS1, y=NMDS2, colour=Distance), linewidth=1) + # adding covariance ellipses according to group
    guides(color = guide_legend(override.aes = list(linetype=c(NA,NA)))) + # removes lines from legend
    theme_bw() + # adding theme
    theme(panel.grid = element_blank(), # remove background grid
          legend.text = element_text(size = 12), # increase legend text size
          legend.title = element_text(size = 14), # increase legend title size
          axis.text = element_text(size = 12), # increase axis text size
          axis.title = element_text(size = 14)) + # increase axis title size) 
    scale_color_manual(name = "Group", # legend title
                       labels = c("area", "flower"), # adjusting legend labels
                       values = c("darkgreen", "orange")) + # customizing colors
    scale_shape_manual("Group", # legend title
                       labels = c("area", "flower"), # adjusting legend labels
                       values = c(17, 19))+ # customizing shapes
    geom_text(aes(label = Site), hjust = 0.5, vjust = -0.5)) # adding site labels

# data analysis
# using a PERMANOVA (PERmutational Multivariate ANalysis Of VAriance) 
# to test the differences in community composition
permanova_group <- adonis2(as.matrix(nmds_data [,3:85]) ~ Group, nmds_data,
                     permutations = 999, method = "bray") 
# permutations is the number of possible arrangements the data could take
# using permutations = 999 is standard practice in science and across the literature
# can increase to get a more thorough analysis
# use method = bray as it is what we previously used to calculate pairwise distances

permanova_group
#p (< 0.01) indicates significant differences in community composition between groups

# check model assumptions
# check for multivariate homogeneity of group variances
# generate distance matrix from pollinator community matrix
dist <- vegdist(community_matrix, method = "bray")

# use betadisper test to check for multivariate homogeneity of group variances
dispersion <- betadisper(dist, group=nmds_data$Group)
permutest(dispersion)
# test has given a non-significant result (p ~ 0.53), meaning that variances
# are homogeneous and model meets the assumptions. Model results can be trusted.


# simper (SIMilarity PERcentage) analysis
# compare community differences, report what species are driving those differences.
simper <- with(nmds_data[,c(1,2)], simper(as.matrix(nmds_data[,3:85]), Group, permutations = 100))

# see most influential species contributing to community differences
simper
#Species represented in this output are responsible for 70% of the observed variation. 
#Each time it adds another taxonomic group it shows its cumulative contribution of 
#rather than its individual one.

summary(simper)


#differences sites----
# data analysis
#PERMANOVA to test differences in community composition between sites
permanova_site <- adonis2(as.matrix(nmds_data [,3:85]) ~ Site, nmds_data,
                     permutations = 999, method = "bray")

permanova_site
#p (~0.25) indicates proportion of permutated F-ratios which are greater or 
#equal to the observed statistic

# check model assumptions
# generate distance matrix from pollinator community matrix
dist <- vegdist(community_matrix, method = "bray")

# use betadisper test to check for multivariate homogeneity of group variances
dispersion_site <- betadisper(dist, group=nmds_data$Site)
permutest(dispersion_site)
# significant result (p < 0.01) -> non-homogeneous variances, model does not 
# meet the assumptions, model results cannot be trusted.

# Include dispersion in model to control for non-homogeneous variance
# Calculate dispersion
dispersion <- betadisper(vegdist(as.matrix(nmds_data[, 3:85]), method = "bray"), nmds_data$Site)

# Extract distances to centroids
dispersion_distances <- dispersion$distances

# Add dispersion distances to data
nmds_data$Dispersion <- dispersion_distances

# Run adonis2 including dispersion as a factor
permanova_site <- adonis2(as.matrix(nmds_data[, 3:85]) ~ Site + Dispersion, data = nmds_data,
                          permutations = 999, method = "bray")

permanova_site
#p (~0.27) indicates no significant community differences between sites
