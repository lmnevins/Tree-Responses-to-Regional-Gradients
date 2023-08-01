######################################################################
#
#  Tree Responses to Regional Gradients -- Mapping Gradients for Figure 1
#
#  L. McKinley Nevins, laura.nevins@wsu.edu, 11 Nov, 2022
#
##      GNU General Public License, Version 3.0    ###################

#Data loading and cleaning for figure 1. There's a lot of excess stuff here but it's
#just easier to include. 

### preamble
require(ecole)
require(data.table)
require(sensitivity)
require(ggplot2)
require(dplyr)
require(tidyverse)
library(maps)
library(mapdata)

### plotting specs
x1 <- expression(log[10] ~ CMD ~ (mm ~ y^{-1}))
x2 <- expression(Fire ~ probability^{1/3})
x3 <- 'Composite stress'
x4 <- 'Funct Neigh x CMD'
x5 <- 'Funct Neigh x Fire'
y2 <- 'Growth (RGR)'
y3 <- 'Survival'
y4 <- 'Crowding (NCI)'
y5 <- 'Funct Neigh (PC1)' # NCIS
y6 <- 'Funct Neigh (PC2)' # NCIS

### load data
setwd("~/Dropbox/Trait_data/Data/Tree_demography_data")
load('pnw_tree-3.rda', verbose=T)
d <- pnw_tree  ;  rm(pnw_tree)
### load functional neighborhood SES per individual
#     NCIS = func neighborhood values =
#        pairwise trait deviations, multiplied by DBH^2, divided by dist^2
#     SES = deviation from random expectation, scaled by SD of random draws
#       problem: should scale and center traits first
#       problem: SES involving congeners always = 0 since traits identical
d[, c('pc1','pc2') := NULL] # rm climate PCs (avoid name conflict w trait PCs)
setwd("~/Dropbox/Trait_data/Results/Functional_neighborhood_SES")
fnm <- c('ses_PC1_alltraits.csv','ses_PC2_alltraits.csv')
ses <- lapply(fnm, fread, colClasses = c('character', 'numeric'))
ses <- Reduce(merge,ses)
names(ses) <- gsub('stand_size_', '', tolower(names(ses)))
names(ses)[names(ses)=='focal'] <- 'tre_cn'
d   <- merge(d, ses, by='tre_cn')
rm(ses)

### exclude juveniles < 5.0 inches (12.7 cm) per GRM and Fenn et al. (2020)
i <- !((!is.na(d$dia_y1) & d$dia_y1 < 5) |
         (!is.na(d$dia_y2) & d$dia_y2 < 5))
d <- d[i,]

### transformations
d$pai     <- d$dia_y2 - d$dia_y1  # periodic annual increment (untransformed)
d$dia_y1  <- log10(1 + d$dia_y1)  # log10-transform
d$dia_y2  <- log10(1 + d$dia_y2)  # log10-transform
d$annincr <- d$dia_y2 - d$dia_y1  # annual increment (transformed)
d$rgr     <- d$annincr / d$dia_y1 # relative growth rate (relative to y1)

### summarize means of SPP in PLOT (SPP/PLOT = sample unit)
means <- d[, list(
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
means <- cbind(means, d[match(means$plt_cn, d$plt_cn), ..j]) ; rm(j) # plot match

p <- means

### trim extreme NCI values > 99.9 percentile
p$plt_nci[p$plt_nci > 2] <- NA
p$rgr[p$rgr > 0.10]      <- NA
is.na(p) <- is.na(p)


#Relationship between cmd and fire? 
hist(p$cmdlog)

hist(p$firelog)

plot(cmdlog ~ firelog, data = p)

#get correlation coefficient
cor(p$cmdlog, p$firelog)

test <- cor.test(p$cmdlog, p$firelog)
test

model <- lm(cmdlog ~ firelog, data = p)

summary(model)

###############################################################
##Mapping the environmental gradients in figure 1

#just cmd 
cc <- ecole::colvec(p$cmdlog, begin=0.1, end=0.95, alpha=0.80)

setwd("~/Dropbox/WSU/Cascade_Project/J of Ecol Submission/Figures")
png(filename='./fig_1A_cmd.png', wid=3,
    hei=4, uni='in', bg='transparent', res=700)

maps::map(database = 'usa', boundary = TRUE, xlim=c(-125,-116.5), ylim=c(38.5,49.5), col = "white", interior = TRUE, fill=FALSE)
points(p$lon, p$lat, pch=19, col=cc, cex=0.5)  #plot my sample sites
map.axes(cex.axis=0.8)


ggplot(p, aes(x = lon, y = lat)) +
  geom_point(color = cc) +
  xlab('Latitude') + 
  ylab('Longitude') +
  theme_classic()

dev.off()


#just fire probability 
cc <- ecole::colvec(p$firelog, begin=0.1, end=0.95, alpha=0.80)

png(filename='./fig_1B_fire.png', wid=3,
    hei=4, uni='in', bg='transparent', res=700)

ggplot(p, aes(x = lon, y = lat)) +
  geom_point(color = cc) +
  xlab('Latitude') + 
  ylab('Longitude') +
  theme_classic()

dev.off()
