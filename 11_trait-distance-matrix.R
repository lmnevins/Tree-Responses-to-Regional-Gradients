######################################################################
#
#  Tree Responses to Regional Gradients -- generating trait distance matrices
#
#  Jenny Zambrano, jenny.zambrano@wsu.edu
#
##      GNU General Public License, Version 3.0    ###################

require(data.table)
require(dplyr)
require(stringr)
require(R.utils)
require(tidyr)
require(tibble)

##################################SINGLE TRAITS#######################################
#Load data 
setwd("~/Dropbox/Cascade_project/Data/Nick_output")
traits_phylo <- read.csv("phy_imp_output.csv")

setwd("~/Dropbox/Cascade_project/Results/PCA results")
fire <- read.csv("PCA_phylopars_scores_fire.csv")

drought <- read.csv("PCA_phylopars_scores_drought.csv")

all_traits <- read.csv("PCA_phylopars_scores_alltraits.csv")

# Add rownames and combine datasets for fire and drought
rownames(traits_phylo) <- traits_phylo$Species

rownames(fire) <- fire$Species
fire <- select(fire, PC1_fire = PC1, PC2_fire = PC2)

rownames(drought) <- drought$Species
drought <- select(drought, PC1_drought = PC1, PC2_drought = PC2)

traits <- cbind(traits_phylo, fire, drought)

# Add rownames for all_traits
rownames(all_traits) <- all_traits$Species
all_traits <- select(all_traits, PC1, PC2)

# Similarity (1 - absolute distance) for fire and drought
sim_traits_abs <- lapply(2:13, function(i) as.matrix(dist(traits[, i, drop = FALSE])))
names(sim_traits_abs) <- colnames(traits)[2:13]

# Similarity (1 - absolute distance) for all traits
sim_all_traits_abs <- lapply(1:2, function(i) as.matrix(dist(all_traits[, i, drop = FALSE])))
names(sim_all_traits_abs) <- colnames(all_traits)[1:2]

# Hierarchical distance for fire and drought
dist_all_traits_hier <- lapply(1:2, function(i) outer(all_traits[, i], all_traits[, i], FUN = "-"))
names(dist_all_traits_hier) <- colnames(all_traits)[1:2]
for (i in 1) {
  rownames(dist_all_traits_hier[[i]]) <- traits_phylo$Species
  colnames(dist_all_traits_hier[[i]]) <- traits_phylo$Species
}

# Save results
setwd("~/Dropbox/Cascade_project/Data/Trait_Data")
saveRDS(sim_traits_abs, "sim_traits_abs.RData")
saveRDS(sim_all_traits_abs, "sim_all_traits_abs.RData")
saveRDS(dist_all_traits_hier, "dist_all_traits_hier.RData")

