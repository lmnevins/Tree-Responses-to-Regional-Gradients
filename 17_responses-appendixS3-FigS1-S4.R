################################################################################################################################
#
#   Tree Responses to Regional Gradients -- Species-specific demographic responses across neighborhoods - Appendix S3 Fig. S1-S4 
#
#   L. McKinley Nevins, laura.nevins@wsu.edu, 4 July 2023
#   Jenny Zambrano, jenny.zambrano@wsu.edu 
#   Rob Smith, phytomosaic@gmail.com
#   
##      GNU General Public License, Version 3.0    #############################################################################

require(ecole)     
require(data.table)
require(dplyr)
require(raster)
require(ggplot2)
require(tidyverse)
require(car)
require(sjPlot)
require(gratia)
require(performance)

##################################################################################

##load in raw data and make survival and growth datasets that don't have any variables scaled

### load data
setwd("~/Documents/McKinley/RStudio/Data")

load('pnw_tree-3.rda', verbose=T)
d         <- pnw_tree  ;  rm(pnw_tree)

### read functional neighborhood SES per individual
d[, c('pc1', 'pc2') := NULL] # rm climate PCs (avoid name conflict w trait PCs)

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


# reassign
trees <- d

# change name format to remove underscore and capitalize genus
trees$spp <- gsub('_', ' ', trees$spp)
trees$spp <- gsub('^(\\w)(\\w+)', '\\U\\1\\L\\2', trees$spp, perl = TRUE)

trees         <- transform(trees, spp = as.factor(spp)) #species as a factor 


### summarize means of SPP in PLOT (SPP/PLOT = sample unit)
p <- trees[, list(
  surv         = mean(surv, na.rm=T),
  rgr          = mean(rgr, na.rm=T),
  plt_nci      = median(i_nci, na.rm=T),
  plt_ses_pc1  = mean(pc1, na.rm=T),
  plt_ses_pc2  = mean(pc2, na.rm=T)
), by=list(plt_cn, spp)] ### <-- (SPP/PLOT = sample unit)
# remove all the INDIVIDUAL-level measures
j <- !names(d) %in% c('tre_cn','dist','subp','subpid','x','y',
                      'is_seed','is_edge','dia_y1','dia_y2','annincr',
                      'plt_cn','spp','qmd','surv','rgr','plt_nci',
                      'i_nci','i_ncilog','pc1','pc2','pai')
p <- cbind(p, d[match(p$plt_cn, d$plt_cn), ..j]) ; rm(j) # plot match

### trim extreme NCI values > 99.9 percentile
p$plt_nci[p$plt_nci > 2] <- NA
p$rgr[p$rgr > 0.10]      <- NA
is.na(p) <- is.na(p)  # force all NaN to NA


# trim down to variables of interest
plot_level <- as.data.frame(dplyr::select(.data=p, plt_cn, lat, lon, spp, cmdlog, firelog, plt_ses_pc1, surv, rgr))

#summarizes for 37,739 plots with data 

# split datasets for growth and survival before removing NA's so data stay intact
survival <- dplyr::select(.data=plot_level, plt_cn, lat, lon, spp, cmdlog, firelog, plt_ses_pc1, surv)
growth   <- dplyr::select(.data=plot_level, plt_cn, lat, lon, spp, cmdlog, firelog, plt_ses_pc1, rgr)

survival <-  tidyr::drop_na(survival) #drops to 35,651 observations = 2,088 NA's
growth   <-  tidyr::drop_na(growth) #drops to 33,886 observations =3,853 NA's

####subset these for our 15 common species of interest 
surv_15 = survival[survival$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                       "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                       "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 

#reduces to 26,683 observations for the 15 most abundant species 

growth_15 = growth[growth$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                     "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                     "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 

#reduces to 25,625 observations for the 15 most abundant species 


####split neighborhoods in both survival and growth datasets 
#when defining the breaks for the categories, if you bump the last number in the decimals for the max and min values, 
#it makes them inclusive and avoids introducing any NAs into the data 

#survival
hist(surv_15$plt_ses_pc1) #distribution skewed pretty far to the left 

mean(surv_15$plt_ses_pc1)
# mean = -0.0001537152

range(surv_15$plt_ses_pc1)
# -1.861373  3.903647

quantile(surv_15$plt_ses_pc1)

#25% = -0.9715723
#75% = 0.9570858


summary(surv_15$plt_ses_pc1)

surv_15$pc1cat <- cut(surv_15$plt_ses_pc1,
                          breaks=c(3.903648, 0.957086, -0.971572, -1.861374), # max, 75%, 25%, min
                          labels=c('Less Dissimilar', 'Intermediate Dissimilarity', 'More Dissimilar'))

summary(surv_15$pc1cat) #most observations are in the intermediate dissimilarity
# less - n=6,671; intermediate - n=13,341; more - n=6,671

#growth 
hist(growth_15$plt_ses_pc1) #distribution skewed pretty far to the left 

mean(growth_15$plt_ses_pc1)
# mean = -0.002256544

range(growth_15$plt_ses_pc1)
# -1.861373  3.903647

quantile(growth_15$plt_ses_pc1)

#25% = -0.9730842
#75% = 0.9538749

summary(growth_15$plt_ses_pc1)


growth_15$pc1cat <- cut(growth_15$plt_ses_pc1,
                            breaks=c(3.903648, 0.953875, -0.973084, -1.861373),
                            labels=c('Less Dissimilar', 'Intermediate Dissimilarity', 'More Dissimilar'))

summary(growth_15$pc1cat) #most observations are in the intermediate dissimilarity
# less - n=6,407; intermediate - n=12,812; more - n=6,406 

#####################################################################################################

#reorder for plotting - this is flipped from the standard order, but necessary to get it to work correctly 
growth_15$spp <- factor(growth_15$spp, levels=rev(c("Juniperus occidentalis", "Calocedrus decurrens", "Quercus chrysolepis", "Abies lasiocarpa", "Alnus rubra", "Thuja plicata", "Tsuga mertensiana", 
                                                            "Lithocarpus densiflorus", "Abies grandis", "Abies amabilis", "Abies concolor", "Pinus contorta", "Tsuga heterophylla", "Pinus ponderosa", "Pseudotsuga menziesii")))


surv_15$spp <- factor(surv_15$spp, levels=rev(c("Juniperus occidentalis", "Calocedrus decurrens", "Quercus chrysolepis", "Abies lasiocarpa", "Alnus rubra", "Thuja plicata", "Tsuga mertensiana", 
                                                    "Lithocarpus densiflorus", "Abies grandis", "Abies amabilis", "Abies concolor", "Pinus contorta", "Tsuga heterophylla", "Pinus ponderosa", "Pseudotsuga menziesii")))


##Figures for Appendix S3 to visualize variation in species-specific survival and growth across the levels of 
#functional neighborhood dissimilarity
##this is using the raw data, not the bootstrapped model results, because they need to be mapped over the real enviro and neighborhood data 

######### Generating Plots########3

##survival
#CMD
surv_cmd_plot <- ggplot(data = surv_15, aes(x = cmdlog, y = surv, color = pc1cat)) +
  geom_point(alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'gray', linewidth = 0.5) +
  facet_wrap( ~ spp) +
  labs(x = "log10 CMD", y = "Survival") +
  theme_bw()

surv_cmd_plot + theme(axis.text.y = element_text(size = 9)) + theme(axis.text.x = element_text(size = 9)) + 
  theme(strip.text = element_text(face = "italic", size = 8)) +
  labs(color = "Functional Neighborhood") + scale_color_manual(values=c("darkorchid4", "maroon", "darkorange"))


#fire
surv_fire_plot <- ggplot(data = surv_15, aes(x = firelog, y = surv, color = pc1cat)) +
  geom_point(alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'gray', linewidth = 0.5) +
  facet_wrap( ~ spp) +
  labs(x = "Wildfire Probability", y = "Survival") +
  theme_bw()

surv_fire_plot + theme(axis.text.y = element_text(size = 9)) + theme(axis.text.x = element_text(size = 9)) + 
  theme(strip.text = element_text(face = "italic", size = 8)) +
  labs(color = "Functional Neighborhood") + scale_color_manual(values=c("darkorchid4", "maroon", "darkorange"))


##growth
#CMD
growth_cmd_plot <- ggplot(data = growth_15, aes(x = cmdlog, y = rgr, color = pc1cat)) +
  geom_point(alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'gray', linewidth = 0.5) +
  facet_wrap( ~ spp) +
  labs(x = "log10 CMD", y = "Relative Growth Rate") +
  theme_bw()

growth_cmd_plot + theme(axis.text.y = element_text(size = 9)) + theme(axis.text.x = element_text(size = 9)) + 
  theme(strip.text = element_text(face = "italic", size = 8)) +
  labs(color = "Functional Neighborhood") + scale_color_manual(values=c("darkorchid4", "maroon", "darkorange"))

#fire
growth_fire_plot <- ggplot(data = growth_15, aes(x = firelog, y = rgr, color = pc1cat)) +
  geom_point(alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'gray', linewidth = 0.5) +
  facet_wrap( ~ spp) +
  labs(x = "Wildfire Probability", y = "Relative Growth Rate") +
  theme_bw()

growth_fire_plot + theme(axis.text.y = element_text(size = 9)) + theme(axis.text.x = element_text(size = 9)) + 
  theme(strip.text = element_text(face = "italic", size = 8)) +
  labs(color = "Functional Neighborhood") + scale_color_manual(values=c("darkorchid4", "maroon", "darkorange"))
