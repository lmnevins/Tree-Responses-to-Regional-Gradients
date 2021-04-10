require(data.table)
require(dplyr)
require(stringr)
require(R.utils)
require(tidyr)
require(tibble)

#Load data
setwd("~/Dropbox/Cascade_project/Results/PCA results")
lambda_hier <- readRDS("dist_all_traits_hier.RData")

setwd("~/Dropbox/Cascade_project/Data/Tree_demography_data")
neighbors <- read.csv("neighbors_FIA.csv")
neighbors$focal <- as.character((neighbors$focal))

# Load trait data
setwd("~/Dropbox/Cascade_project/Data/Nick_output")

traits <- read.csv("phy_imp_output.csv")
traits_dist <- dist(traits[, -1])
traits <- as.matrix(traits_dist)
dimnames(traits) <- dimnames(lambda_hier[[1]])

# Subset species to those with functional traits
neighbors <- filter(neighbors, fspecies %in% rownames(traits),
                    nspecies %in% rownames(traits)) 
neighbors <- arrange(neighbors, fspecies)
neighbors$fspecies <- as.character(neighbors$fspecies)
neighbors$nspecies <- as.character(neighbors$nspecies)

# Functions for neighorhood effects and functional neighborhood for given trait 
trait_hier <- function(neighbors, trait) {
  sum(lambda_hier[[trait]][neighbors$fspecies, neighbors$nspecies] * neighbors$size_dist2)
}

# Function to calculate functional similarity for all traits
functional_all <- function(neighbors) {
  result <- sapply(names(lambda_hier), function (trait) {
    c(trait_hier(neighbors, trait))
  })
  names(result) <- c(paste0("ncih_", names(lambda_hier)))
  data.frame(t(result))
}

#Randomize and save results
setwd("~/Dropbox/Cascade_project/Results/Trait_hierarchy")

neigh_funct <- group_by(neighbors, focal, year) %>% do(functional_all(.))

#Scaling of data
neigh_funct_st <- as.data.frame(scale(log1p(neigh_funct[, 3:14])))

neigh_funct_st[is.na(neigh_funct_st)] <- 0

neigh_all <- bind_cols(neigh_funct[, 1:2], neigh_funct_st)

write.csv(neigh_all, "Trait_hierarchy. csv", row.names = FALSE)
