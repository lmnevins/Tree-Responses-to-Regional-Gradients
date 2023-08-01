################################################################################################################################
#
#   Tree Responses to Regional Gradients -- generate survival and growth datasets for modeling 
#
#   L. McKinley Nevins, laura.nevins@wsu.edu, 16 April 2023
#   Jenny Zambrano, jenny.zambrano@wsu.edu 
#   Rob Smith, phytomosaic@gmail.com
#   
##      GNU General Public License, Version 3.0    #############################################################################


require(ecole)     
require(data.table)
require(dplyr)
require(tidyverse)


### load data
setwd("~/Dropbox/Trait_data/Data/Tree_demography_data")
load('pnw_tree-3.rda', verbose=T)
d         <- pnw_tree  ;  rm(pnw_tree)

### read functional neighborhood SES per individual
d[, c('pc1', 'pc2') := NULL] # rm climate PCs (avoid name conflict w trait PCs)

setwd("~/Dropbox/Trait_data/Results/Functional_neighborhood_SES")

fnm <- c('ses_PC1_alltraits.csv', 'ses_PC2_alltraits.csv')
ses <- lapply(fnm, fread, colClasses = c('character', 'numeric'))
ses <- Reduce(merge,ses)
names(ses) <- gsub('stand_size_', '', tolower(names(ses)))
names(ses)[names(ses)=='focal'] <- 'tre_cn'
d   <- merge(d, ses, by='tre_cn')
rm(ses)

### exclude juveniles < 5.0 inches (12.7 cm) per Fenn et al. (2020)
i <- !((!is.na(d$dia_y1) & d$dia_y1 < 5) |
         (!is.na(d$dia_y2) & d$dia_y2 < 5))
d <- d[i,]
rm(i)

#Goes from 384,089 individuals to 336,428 individuals after juveniles are excluded 

# transformations
d$pai     <- d$dia_y2 - d$dia_y1  # annual increment (untransformed)
d$dia_y1  <- log10(1 + d$dia_y1)  # log10-transform
d$dia_y2  <- log10(1 + d$dia_y2)  # log10-transform
d$annincr <- d$dia_y2 - d$dia_y1  # annual increment (transformed)
d$rgr     <- d$annincr / d$dia_y1 # relative growth rate (relative to y1)
d         <- transform(d, censusyr = as.integer(measyr_2)) #census year as an integer 
d$census_int   <- d$measyr_2 - d$measyr_1 #calculate a sampling interval that takes into account different time spans between census years 1 and 2
d   <- transform(d, census_int = as.integer(census_int)) #census year as an integer 

# scaling
#Including scaling of the diameter at year 1 for each tree, to get a standardized DBH
j <- c('cmdlog','firelog','pc1','dia_y1') # ,'surv','rgr' - cut out if we don't want to scale responses
d[, (j) := lapply(.SD, scale), .SDcols=j]

# reassign
trees <- d

# change name format to remove underscore and capitalize genus
trees$spp <- gsub('_', ' ', trees$spp)
trees$spp <- gsub('^(\\w)(\\w+)', '\\U\\1\\L\\2', trees$spp, perl = TRUE)

trees         <- transform(trees, spp = as.factor(spp)) #species as a factor 

# trim down to variables of interest
trees <- as.data.frame(dplyr::select(.data=trees, tre_cn, plt_cn, lat, lon, spp, cmdlog, firelog, pc1,
                                     censusyr, dia_y1, surv, rgr, census_int))

# split datasets for growth and survival before removing NA's so data stay intact
survival <- dplyr::select(.data=trees, tre_cn, plt_cn, lat, lon, spp, cmdlog, firelog, pc1,
                          censusyr, dia_y1, surv, census_int)
growth   <- dplyr::select(.data=trees, tre_cn, plt_cn, lat, lon, spp, cmdlog, firelog, pc1,
                          censusyr, dia_y1, rgr, census_int)
survival <-  tidyr::drop_na(survival) #29,407 rows removed due to NA's in PC1 and PC2 data (58,701 total observations missing)
growth   <-  tidyr::drop_na(growth) # 65,949 rows removed due to NA's in PC1 and PC2 data, and zeros in rgr from indivs with 0 survival (98,200 total observations missing)


#save datafiles to be used on Kamiak 
setwd("~/Dropbox/WSU/Cascade_Project/Kamiak/Resub2")
write.csv(survival, "survival_full.csv", row.names=FALSE)
write.csv(growth, "growth_full.csv", row.names=FALSE)
