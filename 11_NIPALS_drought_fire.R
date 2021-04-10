# NIPALS
# Using drought and fire trait data from FIA and TRY
# Based on script from Rob Smith
#
# McKinley Nevins, 21 June 2020
#
# Explanation of what the NIPALS function does from https://stats.stackexchange.com/questions/35561/?nipals
# "the NIPALS algorithm interpolates the missing point using a least squares fit but give 
# the missing data no influence on the model. Successive iterations refine the missing value
# by simply multiplying the score and the loading for that point." 
# Stands for Non-linear Iterative Partial Least Squares

install.packages("ps")
library(githubinstall)
githubinstall("ecole")
require(data.table)
require(ecole)
require(ps)
require(vegan)
require(viridis)
require(dplyr)

#Package with the nipals function
require(nipals)


setwd("~/Dropbox/Trait_Data//Results/PCA results")
### load USDA and FIA traits
trf <- read.csv("Fixed_all_fire_traits.csv", stringsAsFactors=F)
setwd("~/Dropbox/Trait_Data/Data/Trait_Data")
trd <- read.csv("Species_drought_traits.csv", stringsAsFactors=F)
# The fire dataset has three more species in it than the drought dataset, because there are only 53 species 
# which have data for at least one of the drought traits

# We are going to consider the drought and fire traits separately
# These actions change the way the columns are titles and the species column is displayed
#Fire
trf$Species <- ecole::clean_text(trf$Species, T)
names(trf)  <- ecole::clean_text(names(trf), T)
dimnames(trf)[[1]] <- trf$species
trf$species <- NULL
trf$seed_mass[trf$seed_mass > 500000] <- NA # fix one crazy value

#Drought
trd$Species <- ecole::clean_text(trd$Species, T)
names(trd)  <- ecole::clean_text(names(trd), T)
dimnames(trd)[[1]] <- trd$species
trd$species <- NULL

#-----------------------------------------------------------------
# NIPALS PCA for fire
(pcfire <- ade4::nipals(trf, nf=3, rec=T, niter=199))
tras_fire <- pcfire$rec # keep the *imputed* z-score values that have filled in the NAs
scr_fire  <- pcfire$li  # keep the species scores
pcfire$eig         # eigenvalues
dimnames(scr_fire)[[2]] <- paste0('PC',1:3) # adds species names to the species scores

# Save scores
write.csv(scr_fire, "Fire_traits_scores_nipals.csv", row.names = TRUE)

### labels to overlay hardwood/conifer
is_coniffire <- rep(TRUE,NROW(trf))
is_coniffire[c(7:13,17,18,20,24,25,38:41,43:49,56)] <- FALSE

### PCA biplot
plot(scr_fire*1.25, type='n', main="Fire Traits", ylab='PC2', xlab='PC1',
     xlim=c(-7, 20), ylim=c(-5,4))
points(scr_fire, pch=3, col=c(1,2)[(is_coniffire) + 1])
text(scr_fire, lab=dimnames(scr_fire)[[1]], cex=0.5, pos=4,
     col=c(1,2)[(is_coniffire) + 1])
arrows(0,0,pcfire$co[,1]*5,pcfire$co[,2]*5,col='#0000FF90',
       lwd=2,ang=15,len=0.2)
text(pcfire$co*5, names(trf), col='#0000FFDD', cex=1, pos=4, font=2)
legend('bottomleft',leg=c('Hardwood','Conifer'),
       fill=c(1,2),bty='n',cex=0.7)


#NIPALS PCA for drought 
#Had to increase the number of iterations, because it was reaching the maximum before 
#completing the first axis
(pcdrought <- ade4::nipals(trd, nf=3, rec=T, niter=300))
tras_drought <- pcdrought$rec # keep the *imputed* z-score values that have filled in the NAs
scr_drought  <- pcdrought$li  # keep the species scores
pcdrought$eig         # eigenvalues
dimnames(scr_drought)[[2]] <- paste0('PC',1:3) # adds species names to the species scores

# Save scores
write.csv(scr_drought, "Drought_traits_scores_nipals.csv", row.names = TRUE)

### labels to overlay hardwood/conifer
is_conifdrought <- rep(TRUE,NROW(trd))
is_conifdrought[c(7:13,17,18,19,23,24,36:39,41:47,53)] <- FALSE

### PCA biplot
plot(scr_drought*1.25, type='n', main="Drought Traits", ylab='PC2', xlab='PC1',
     xlim=c(-5, 7), ylim=c(-3,3))
points(scr_drought, pch=3, col=c(1,2)[(is_conifdrought) + 1])
text(scr_drought, lab=dimnames(scr_drought)[[1]], cex=0.5, pos=4,
     col=c(1,2)[(is_conifdrought) + 1])
arrows(0,0,pcdrought$co[,1]*5,pcdrought$co[,2]*5,col='#0000FF90',
       lwd=2,ang=15,len=0.2)
text(pcdrought$co*5, names(trd), col='#0000FFDD', cex=1, pos=4, font=2)
legend('bottomleft',leg=c('Hardwood','Conifer'),
       fill=c(1,2),bty='n',cex=0.7)

#---------------------------------------------------------
# NIPALS PCA for Fire Excluding SLA

(pcfirenoSLA <- ade4::nipals(trf[2:5], nf=3, rec=T, niter=199))
tras_firenoSLA <- pcfirenoSLA$rec # keep the *imputed* z-score values that have filled in the NAs
scr_firenoSLA  <- pcfirenoSLA$li  # keep the species scores
pcfirenoSLA$eig         # eigenvalues
dimnames(scr_firenoSLA)[[2]] <- paste0('PC',1:3) # adds species names to the species scores

# Save scores
setwd("~/Dropbox/Trait_Data/Results/PCA Results")
write.csv(scr_firenoSLA, "Fire_traits_scores_nipals_noSLA.csv", row.names = TRUE)

### labels to overlay hardwood/conifer
is_coniffire <- rep(TRUE,NROW(trf))
is_coniffire[c(7:13,17,18,20,24,25,38:41,43:49,56)] <- FALSE

### PCA biplot
plot(scr_firenoSLA*1.25, type='n', main="Fire Traits Excluding SLA", ylab='PC2', xlab='PC1',
     xlim=c(-20, 15), ylim=c(-5,4))
points(scr_firenoSLA, pch=3, col=c(1,2)[(is_coniffire) + 1])
text(scr_firenoSLA, lab=dimnames(scr_firenoSLA)[[1]], cex=0.5, pos=4,
     col=c(1,2)[(is_coniffire) + 1])
arrows(0,0,pcfirenoSLA$co[,1]*5,pcfirenoSLA$co[,2]*5,col='#0000FF90',
       lwd=2,ang=15,len=0.2)
text(pcfirenoSLA$co*5, names(trf[2:5]), col='#0000FFDD', cex=1, pos=4, font=2)
legend('bottomleft',leg=c('Hardwood','Conifer'),
       fill=c(1,2),bty='n',cex=0.7)

#NIPALS PCA for drought excluding SLA
(pcdroughtnoSLA <- ade4::nipals(trd[1:3], nf=3, rec=T, niter=800))
tras_droughtnoSLa <- pcdroughtnoSLA$rec # keep the *imputed* z-score values that have filled in the NAs
scr_droughtnoSLA  <- pcdroughtnoSLA$li  # keep the species scores
pcdroughtnoSLA$eig         # eigenvalues
dimnames(scr_droughtnoSLA)[[2]] <- paste0('PC',1:3) # adds species names to the species scores

# Save scores
write.csv(scr_droughtnoSLA, "Drought_traits_scores_nipals_noSLA.csv", row.names = TRUE)

### labels to overlay hardwood/conifer
is_conifdrought <- rep(TRUE,NROW(trd))
is_conifdrought[c(7:13,17,18,19,23,24,36:39,41:47,53)] <- FALSE

### PCA biplot
plot(scr_droughtnoSLA*1.25, type='n', main="Drought Traits Excluding SLA", ylab='PC2', xlab='PC1',
     xlim=c(-5, 7), ylim=c(-3,4))
points(scr_droughtnoSLA, pch=3, col=c(1,2)[(is_conifdrought) + 1])
text(scr_droughtnoSLA, lab=dimnames(scr_droughtnoSLA)[[1]], cex=0.5, pos=4,
     col=c(1,2)[(is_conifdrought) + 1])
arrows(0,0,pcdroughtnoSLA$co[,1]*5,pcdroughtnoSLA$co[,2]*5,col='#0000FF90',
       lwd=2,ang=15,len=0.2)
text(pcdroughtnoSLA$co*5, names(trd[1:3]), col='#0000FFDD', cex=1, pos=4, font=2)
legend('bottomleft',leg=c('Hardwood','Conifer'),
       fill=c(1,2),bty='n',cex=0.7)


########################################################################################

# Bootstrap confidence intervals calculations
# Adapted from Rob's script

###Fire Traits First########

### bootstrap NIPALS
B <- 5000 # number of bootstraps (increase >999 for precision)
`ffire` <- function(x, ...) {
  i   <- sample(1:NROW(x), size=NROW(x), replace=T)
  pcfireboot  <- nipals::nipals(trf, ncomp=2, fitted=T,
                        maxit=199, tol=1e-3) # liberal...
  scrfireboot <- as.matrix(cbind(pcfireboot$scores, pcfireboot$fitted))
  dimnames(scrfireboot)[[1]] <- gsub('\\..*', '', dimnames(scrfireboot)[[1]])
  scrfireboot <- scrfireboot[!duplicated(scrfireboot),]
  scrfireboot <- scrfireboot[order(dimnames(scrfireboot)[[1]]),]
  return(scrfireboot)
}
bfire <- lapply(1:B, function(i) {  # ! ! ! TIMEWARN ! ! !
  cat(paste0(i, ifelse(i %% 25 == 0,'\n',' ')))
  ffire(trf)
})
bbfire <- do.call('rbind', bfire)

### 95% CIs for each species 
fireconf <- do.call(data.frame, aggregate(bbfire, list(spp=rownames(bbfire)),
          function(i) quantile(i,prob=c(0.025,0.975))))

### Mean of imputed trait values across all of the repetitions
fireboot <- do.call(data.frame, aggregate(bbfire, list(spp=rownames(bbfire)),
                                          function(i) mean(i)))

## Save results
setwd("~/Dropbox/Trait_Data/Results/PCA Results")
write.csv(fireconf, "Nipals_boot_fire_CI.csv", row.names = TRUE)
write.csv(fireboot, "Nipals_boot_fire.csv", row.names = TRUE) 


###Fire excluding SLA##################

#When SLA is excluded from the fire traits file, one species is all NAs and must be excluded
trfSLA <- select(trf, bark_thickness, litter_decomposition, resprouting, seed_mass)

rownames(trfSLA)[39]->remove
trfSLA<- trfSLA[setdiff(rownames(trfSLA), remove), ]

### bootstrap NIPALS
B <- 5000 # number of bootstraps (increase >999 for precision)
`ffireSLA` <- function(x, ...) {
  i2   <- sample(1:NROW(x), size=NROW(x), replace=T)
  pcfirebootSLA  <- nipals::nipals(trfSLA, ncomp=2, fitted=T,
                                maxit=199, tol=1e-3) # liberal...
  scrfirebootSLA <- as.matrix(cbind(pcfirebootSLA$scores, pcfirebootSLA$fitted))
  dimnames(scrfirebootSLA)[[1]] <- gsub('\\..*', '', dimnames(scrfirebootSLA)[[1]])
  scrfirebootSLA <- scrfirebootSLA[!duplicated(scrfirebootSLA),]
  scrfirebootSLA <- scrfirebootSLA[order(dimnames(scrfirebootSLA)[[1]]),]
  return(scrfirebootSLA)
}
bfireSLA <- lapply(1:B, function(i2) {  # ! ! ! TIMEWARN ! ! !
  cat(paste0(i2, ifelse(i2 %% 25 == 0,'\n',' ')))
  ffireSLA(trfSLA)
})
bbfireSLA <- do.call('rbind', bfireSLA)

### 95% CIs for each species 
fireconfSLA <- do.call(data.frame, aggregate(bbfireSLA, list(spp=rownames(bbfireSLA)),
                                          function(i2) quantile(i2,prob=c(0.025,0.975))))

### Mean of imputed trait values across all of the repetitions
firebootSLA <- do.call(data.frame, aggregate(bbfireSLA, list(spp=rownames(bbfireSLA)),
                                          function(i2) mean(i2)))

## Save results
setwd("~/Dropbox/Trait_Data/Results/PCA Results")
write.csv(fireconfSLA, "Nipals_boot_firenoSLA_CI.csv", row.names = TRUE)
write.csv(firebootSLA, "Nipals_boot_firenoSLA.csv", row.names = TRUE) 


##Drought traits######################################

### bootstrap NIPALS
B <- 5000 # number of bootstraps (increase >999 for precision)
`fdrought` <- function(x, ...) {
  i3   <- sample(1:NROW(x), size=NROW(x), replace=T)
  pcdroughtboot  <- nipals::nipals(trd, ncomp=2, fitted=T,
                                maxit=199, tol=1e-3, startcol = 0) # liberal...
  scrdroughtboot <- as.matrix(cbind(pcdroughtboot$scores, pcdroughtboot$fitted))
  dimnames(scrdroughtboot)[[1]] <- gsub('\\..*', '', dimnames(scrdroughtboot)[[1]])
  scrdroughtboot <- scrdroughtboot[!duplicated(scrdroughtboot),]
  scrdroughtboot <- scrdroughtboot[order(dimnames(scrdroughtboot)[[1]]),]
  return(scrdroughtboot)
}
bdrought <- lapply(1:B, function(i3) {  # ! ! ! TIMEWARN ! ! !
  cat(paste0(i3, ifelse(i3 %% 25 == 0,'\n',' ')))
  fdrought(trd)
})

bbdrought <- do.call('rbind', bdrought)

### 95% CIs for each species 
droughtconf <- do.call(data.frame, aggregate(bbdrought, list(spp=rownames(bbdrought)),
                                          function(i3) quantile(i3,prob=c(0.025,0.975))))

### Mean of imputed trait values across all of the repetitions
droughtboot <- do.call(data.frame, aggregate(bbdrought, list(spp=rownames(bbdrought)),
                                          function(i3) mean(i3)))

## Save results
setwd("~/Dropbox/Trait_Data/Results/PCA Results")
write.csv(droughtconf, "Nipals_boot_drought_CI.csv", row.names = TRUE)
write.csv(droughtboot, "Nipals_boot_drought.csv", row.names = TRUE) 


###Drought no SLA###########################

### bootstrap NIPALS
B <- 5000 # number of bootstraps (increase >999 for precision)
`fdroughtSLA` <- function(x, ...) {
  i4   <- sample(1:NROW(x), size=NROW(x), replace=T)
  pcdroughtbootSLA  <- nipals::nipals(trd[1:3], ncomp=2, fitted=T,
                                   maxit=199, tol=1e-3) # liberal...
  scrdroughtbootSLA <- as.matrix(cbind(pcdroughtbootSLA$scores, pcdroughtbootSLA$fitted))
  dimnames(scrdroughtbootSLA)[[1]] <- gsub('\\..*', '', dimnames(scrdroughtbootSLA)[[1]])
  scrdroughtbootSLA <- scrdroughtbootSLA[!duplicated(scrdroughtbootSLA),]
  scrdroughtbootSLA <- scrdroughtbootSLA[order(dimnames(scrdroughtbootSLA)[[1]]),]
  return(scrdroughtbootSLA)
}
bdroughtSLA <- lapply(1:B, function(i4) {  # ! ! ! TIMEWARN ! ! !
  cat(paste0(i4, ifelse(i4 %% 25 == 0,'\n',' ')))
  fdroughtSLA(trd[1:3])
})
bbdroughtSLA <- do.call('rbind', bdroughtSLA)

### 95% CIs for each species 
droughtconfSLA <- do.call(data.frame, aggregate(bbdroughtSLA, list(spp=rownames(bbdroughtSLA)),
                                             function(i4) quantile(i4,prob=c(0.025,0.975))))

### Mean of imputed trait values across all of the repetitions
droughtbootSLA <- do.call(data.frame, aggregate(bbdroughtSLA, list(spp=rownames(bbdroughtSLA)),
                                             function(i4) mean(i4)))

## Save results
setwd("~/Dropbox/Trait_Data/Results/PCA Results")
write.csv(droughtconfSLA, "Nipals_boot_droughtnoSLA_CI.csv", row.names = TRUE)
write.csv(droughtbootSLA, "Nipals_boot_droughtnoSLA.csv", row.names = TRUE) 




####################################################################################
### plot species replicates in PCA trait space
###   (but first remove some outliers > 4 SD)
o1 <- ecole::outlier_uni(bbfire[,1], mult=4)
o2 <- ecole::outlier_uni(bbfire[,2], mult=3.5)
o  <- c(o1$upper,o1$lower,o2$upper,o2$lower)
rfire  <- vegan::scores(bbfire[-o,1:2])
rnfire <- as.factor(rownames(rfire))
table(rnfire)
png('./fig_00.png', wid=4.5, hei=4.5, uni='in',
    bg='transparent', res=500)
vegan::ordiplot(rfire,type='n',display='sites', bty='l',las=1,
                xaxt='none',yaxt='none', ylab='PC2',xlab='PC1')
ccfire <- viridis::inferno(min(nlevels(rnfire), 99), alpha=0.7)
# points(r, pch=16, col=cc[rn])
vegan::ordispider(rfire, groups=rnfire, col=ccfire)
dev.off()

####    END    ####