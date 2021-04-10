# PCA from Rphylopars imputed data
# Drought and Fire Traits
#
# McKinley Nevins 
# July 22, 2020
#
#####################################

# PCAs from the imputed dataset produced using the phylopars function in the Rphylopars package.

setwd("~/Dropbox/Trait_Data/Data/Nick_output")

library(dplyr)
require(tidyr)
library(vegan)
library(rcompanion)
require(ggplot2)
require(ggfortify)
library(oompaBase)
library(oompaData)
library(mclust)
library(ClassDiscovery)
library(kernlab)
library(changepoint)
library(cpm)
library(PCDimension)

# Load imputed data file of all traits 
alltraits <- read.csv("phy_imp_output.csv")


#Split up data set for drought and fire traits 
drought <- select(alltraits, Species, SLA, P50, Stem_conductivity, Root_depth)
fire <- select(alltraits, Species, Bark_thickness, Litter_decomposition, Seed_mass, Resprouting)

####################################
#PCA's
#All Traits
alltraits.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(alltraits[2:9])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$Species = alltraits$Species
summary(alltraits.pca)

# Save loadings for all traits
setwd("~/Dropbox/Trait_Data/Results/PCA results")
write.csv(loadings.alltraits, "PCA_phylopars_loadings_alltraits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "PCA_phylopars_scores_alltraits.csv")

#Biplot
plot <- ggplot2::autoplot(alltraits.pca, data = alltraits, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.5,0.6), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for All Traits")+
  theme_classic()

print(plot)

#Save biplot
ggsave("PCA_phylopars_alltraits.pdf", width = 6, height = 5, units = "in")

##########################################
#Drought Traits
drought.pca = prcomp(drought[2:5], center = T, scale = T)

sd.drought = drought.pca$sdev
loadings.drought = drought.pca$rotation
trait.names.drought = colnames(drought[2:5])
scores.drought = as.data.frame(drought.pca$x)
scores.drought$Species = drought$Species
summary(drought.pca)

# Save loadings for drought traits
write.csv(loadings.drought, "PCA_phylopars_loadings_drought.csv", row.names = TRUE)

#Save species scores
write.csv(scores.drought, "PCA_phylopars_scores_drought.csv")

#######################################
#Fire Traits
fire.pca = prcomp(fire[2:5], center = T, scale = T)

sd.fire = fire.pca$sdev
loadings.fire = fire.pca$rotation
trait.names.fire = colnames(fire[2:5])
scores.fire = as.data.frame(fire.pca$x)
scores.fire$Species = fire$Species
summary(fire.pca)

# Save loadings for fire traits
write.csv(loadings.fire, "PCA_phylopars_loadings_fire.csv", row.names = TRUE)

#Save species scores
write.csv(scores.fire, "PCA_phylopars_scores_fire.csv")




##Broken-Stick test for the significance of the loadings

bs_pca <- alltraits.pca(t(ranData))
