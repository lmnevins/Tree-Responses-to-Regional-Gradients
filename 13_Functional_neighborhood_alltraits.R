require(data.table)
require(dplyr)
require(stringr)
require(R.utils)
require(tidyr)
require(tibble)

#Load data
setwd("~/Dropbox/Trait_Data")
lambda_abs <- readRDS("sim_all_traits_abs.RData")

setwd("~/Dropbox/Cascade_project/Data/Tree_demography_data")
neighbors <- read.csv("neighbors_FIA.csv")
neighbors$focal <- as.character((neighbors$focal))

# Load trait data
setwd("~/Dropbox/Cascade_project/Data/Nick_output")

traits <- read.csv("phy_imp_output.csv")
traits_dist <- dist(traits[, -1])
traits <- as.matrix(traits_dist)
dimnames(traits) <- dimnames(lambda_abs[[1]])

# Subset species to those with functional traits
neighbors <- filter(neighbors, fspecies %in% rownames(traits),
                    nspecies %in% rownames(traits)) 
neighbors <- arrange(neighbors, fspecies)
neighbors$fspecies <- as.character(neighbors$fspecies)
neighbors$nspecies <- as.character(neighbors$nspecies)

# Functions for neighorhood effects and functional neighborhood for given trait 
functional_abs <- function(neighbors, trait) {
  sum(lambda_abs[[trait]][neighbors$fspecies, neighbors$nspecies] * neighbors$size_dist2) 
}


# Function to calculate functional similarity for all traits
functional_all <- function(neighbors) {
  result <- sapply(names(lambda_abs), function (trait) {
    c(functional_abs(neighbors, trait))
  })
  names(result) <- c(paste0("ncis_", names(lambda_abs)))
  data.frame(t(result))
}

#Randomize and save results
setwd("~/Dropbox/Functional_neighborhood_randomizations")

for (i in 1:1000) {
  # Randomize trait differences
  for (g in 1:length(lambda_abs)) {
    rand <- sample(rownames(lambda_abs[[g]]), nrow(lambda_abs[[g]]))
    rownames(lambda_abs[[g]]) <- rand
    colnames(lambda_abs[[g]]) <- rand
  }
  
  neigh_funct <- group_by(neighbors, focal, year) %>% do(functional_all(.))
  
  write.csv(neigh_funct, paste0("functi_neigh", i, ".csv"), row.names = FALSE)    
}
