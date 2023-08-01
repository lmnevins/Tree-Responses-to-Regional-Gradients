######################################################################
#
#  Tree Responses to Regional Gradients -- trait imputation with Rphylopars
#
#
##      GNU General Public License, Version 3.0    ###################


require(Rphylopars)
require(phytools) 
require(data.table)
require(dplyr)
require(stringr)
require(R.utils)
require(tidyr)
require(tibble)

#Load data
setwd("~/Dropbox//Trait_Data/Data/Trait_Data")
traits <- read.csv("trait_data_combine.csv")

setwd("~/Dropbox/Trait_Data/Data/Tree_demography_data")
load(file = "pnw_tree.rda")

#Load phylogeny
setwd("~/Dropbox/Trait_data/Data/phylogeny")
load("pnw_phy.rda")


#Get species names from FIA dataset
pnw_tree$spp <- gsub("_", " ", pnw_tree$spp)
pnw_tree$spp <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
                     pnw_tree$spp, perl = TRUE)

species_fia <- distinct(pnw_tree, spp) %>%
  arrange(spp)

#Keep species from FIA dataset
traits_sub <- arrange(traits, TraitName) %>%
  filter(TraitName != '') %>%
  filter(!is.na(OrigValueStr)) %>%
  filter(SpeciesName %in% species_fia$spp)

traits_unique <- distinct(traits_sub, TraitName)

#Get trait mean
traits_mat<- filter(traits_sub, TraitName %in% c("bark_vol_pct", "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded",
                                                 "resprout_ability","Seed dry mass", "Litter decomposition rate", "Xylem hydraulic vulnerability, xylem cavitation vulnerability, embolism vulnerability, (P20, P50, P80)",
                                                 "Wood (sapwood) specific conductivity (stem specific conductivity)", "root_depth_minimum_inches")) %>%
  
  group_by(SpeciesName, TraitName) %>%
  mutate(OrigValueStr = as.numeric(as.character(OrigValueStr))) %>%
  dplyr::summarize(trait_mean = mean(OrigValueStr))

#Reorganize data to separate each trait
trait_mat_spread <- spread(traits_mat, TraitName, trait_mean, fill = NA, convert = FALSE, drop = TRUE, sep = NULL)
colnames(trait_mat_spread) <- c("species", "Bark_thickness", "SLA", "Litter_decomposition", "Resprouting", "Seed_mass", "P50", "Stem_conductivity", "Root_depth")

trait_mat_spread <- as.data.frame(trait_mat_spread)

#Create trait covariance matrix
trait_mat <- as.matrix(trait_mat_spread[,2:9])
rownames(trait_mat) <- trait_mat_spread$species

trait_mat[is.na(trait_mat)] <- 0

trait_cov <- cov(trait_mat)

#Load Tree
tree <- pnw_phy

#Remove underscore in tip labels and capitalize first character
tree$tip.label <- gsub("_", " ", tree$tip.label)%>%
R.utils::capitalize(tip.label)

#Simulate missing data, I think nmissing should be equal to 0 as we don't want any missing data
sim_data <- simtraits(v = trait_cov,tree = tree,nmissing = 0)
p_BM <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree)

#Original data with missing values
sim_data$trait_data[,2:9]

#View the imputed missing data and imputation variance, from which we construct 95% confidence intervals
p_BM$anc_recon[1:56,] # Data with imputed species means
p_BM$anc_var[1:56,] # Variances for each estimate
p_BM$anc_recon[1:56,] - sqrt(p_BM$anc_var[1:56,])*1.96 # Lower 95% CI
p_BM$anc_recon[1:56,] + sqrt(p_BM$anc_var[1:56,])*1.96 # Upper 95% CI

#Plot phylogeny
plot.phylo(reorder(tree,"postorder"))
nodelabels()

