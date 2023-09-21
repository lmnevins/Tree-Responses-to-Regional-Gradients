####################################################################################
#
#  Tree Responses to Regional Gradients -- Regressions for Figure 3 
#
#    L. McKinley Nevins, laura.nevins@wsu.edu, 24 Nov 2021
#    Rob Smith, phytomosaic@gmail.com
#
#
##      GNU General Public License, Version 3.0    ###################

### preamble
require(ecole)     
require(data.table)
require(dplyr)
require(raster)
require(ggplot2)
require(tidyverse)
require(car)
require(gratia)
require(mgcv)
require(ggpubr)


setwd("~/Documents/McKinley/RStudio/Data")


#######standard data set-up to get dataset of all trees for plotting#########
### load data
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
d         <- transform(d, censusyr = as.integer(measyr_2)) #census year as an integer 
d$census_int   <- d$measyr_2 - d$measyr_1 #calculate a sampling interval that takes into account different time spans between census years 1 and 2
d   <- transform(d, census_int = as.integer(census_int)) #census year as an integer 

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


# Next step is to trim down the data set to individual data for just variables of interest 
trees <- p %>% dplyr::select(plt_cn, lat, lon, spp, cmdlog, firelog, plt_nci, plt_ses_pc1, plt_ses_pc2)

### trim extreme NCI values > 99.9 percentile
trees$plt_nci[trees$plt_nci > 2] <- NA
is.na(trees) <- is.na(trees)  # force all NaN to NA

#convert to dataframe from data table 
trees <- as.data.frame(trees)
class(trees)

#remove na's related to the PCA data 
trees <-  tidyr::drop_na(trees) #removes
#now 35,503 plots 

#################################################################################################################
#Simple Linear Regressions of cmd (cmdlog) and wildfire probability (firelog) with variables for Fig. 3. 

##PC1 functional neighborhood 
#cmdlog --- adjusted r squared 0.184
cmdxpc1_plot <- ggplot(trees, aes(x = cmdlog, y = plt_ses_pc1)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(1.4, 2.5) +
  xlab('log10 CMD') + 
  ylab('PC1 Functional Neighborhood') +
  theme_classic()

cmdxpc1_plot

cmdxpc1_mod <- lm(plt_ses_pc1 ~ cmdlog, data = trees)
summary(cmdxpc1_mod)


#fireprob --- adjusted r squared 0.0868
firexpc1_plot <- ggplot(trees, aes(x = firelog, y = plt_ses_pc1)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(0.0, 0.35) +
  xlab('Wildfire Probability') + 
  ylab('PC1 Functional Neighborhood') +
  theme_classic()

firexpc1_plot

firexpc1_mod <- lm(plt_ses_pc1 ~ firelog, data = trees)
summary(firexpc1_mod)


##PC2 functional neighborhood 
#cmdlog --- adjusted r squared 0.0009
cmdxpc2_plot <- ggplot(trees, aes(x = cmdlog, y = plt_ses_pc2)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(1.4, 2.5) +
  xlab('log10 CMD') + 
  ylab('PC2 Functional Neighborhood') +
  theme_classic()

cmdxpc2_plot

cmdxpc2_mod <- lm(plt_ses_pc2 ~ cmdlog, data = trees)
summary(cmdxpc2_mod)


#fireprob --- adjusted r squared 0.0017
firexpc2_plot <- ggplot(trees, aes(x = firelog, y = plt_ses_pc2)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(0.0, 0.35) +
  xlab('Wildfire Probability') + 
  ylab('PC2 Functional Neighborhood') +
  theme_classic()

firexpc2_plot

firexpc2_mod <- lm(plt_ses_pc2 ~ firelog, data = trees)
summary(firexpc2_mod)


####crowding index

#cmdlog --- adjusted r squared  0.0000
cmdxnci_plot <- ggplot(trees, aes(x = cmdlog, y = plt_nci)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(1.4, 2.5) +
  xlab('log10 CMD') + 
  ylab('Neighborhood Crowding Index') +
  theme_classic()

cmdxnci_plot

cmdxnci_mod <- lm(plt_nci ~ cmdlog, data = trees)
summary(cmdxnci_mod)


#fireprob --- adjusted r squared 0.0000
firexnci_plot <- ggplot(trees, aes(x = firelog, y = plt_nci)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth(method = lm, color = 'red', se = TRUE, fill = 'lightblue') +
  xlim(0.0, 0.35) +
  xlab('Wildfire Probability') + 
  ylab('Neighborhood Crowding Index') +
  theme_classic()

firexnci_plot

firexnci_mod <- lm(plt_nci ~ firelog, data = trees)
summary(firexnci_mod)



# The PC1 neighborhoods are still the only ones that show any 
# significant relationships with the environmental variables in 
# the simple linear regressions 


ggarrange(cmdxpc1_plot, firexpc1_plot, cmdxpc2_plot, firexpc2_plot, cmdxnci_plot, firexnci_plot, nrow=3, ncol=2, common.legend = FALSE, align = "hv", hjust = -1)
