######################################################################
#
#  Tree Responses to Regional Gradients -- determining neighbors
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

#Load data
setwd("~/Dropbox/Cascade_project/Data/Tree_demography_data")
load(file = "pnw_tree.rda")

#Load trait data
#drought <-
#fire <-
#pc_scores_drought <-
#pc_scores_fire <-

#Subset data and calculate relative growth rate
tree_sub <- select(pnw_tree, plt_cn, tre_cn, spp, measyr_1, measyr_2, dia_y1, dia_y2, surv) %>%
            mutate(rgr = (log(dia_y2) - log(dia_y1))/ (measyr_2 - measyr_1))

#Determine focal tree and its neighbors. 
neighbors <- inner_join(tree_sub, tree_sub, by = c("plt_cn","measyr_1")) %>%
  filter(tre_cn.x != tre_cn.y) %>%
  dplyr::select(plot = plt_cn, focal = tre_cn.x, fspecies = spp.x, 
                nspecies = spp.y, nsize = dia_y1.y, year = measyr_1)

#Determine the distance to neighbors
neighbors$distance <- 1

#Remove missing values
neighbors <- filter(neighbors, !is.na(nsize))

# Calculate neighborhood crowding index
neighbors$size_dist2 <- (neighbors$nsize^2 / neighbors$distance^2) 

#Modify species names to match trait data 
neighbors$fspecies <- gsub("_", " ", neighbors$fspecies)
neighbors$fspecies <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
                           neighbors$fspecies, perl = TRUE)

neighbors$nspecies <- gsub("_", " ", neighbors$nspecies)
neighbors$nspecies <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",
                           neighbors$nspecies, perl = TRUE)

write.csv(neighbors, "neighbors_FIA.csv", row.names = FALSE)

