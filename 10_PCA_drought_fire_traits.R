# PCAs for Drought and Fire Traits
#
#
#
# McKinley Nevins 
#
#
# An error was noted in the fire trait data files, so everything was redone with corrected fire trait data
# Fire data set names were kept the same to reduce changes necessary to the script

setwd("~/Dropbox/Trait_Data/Data/Trait_Data")

library(dplyr)
require(tidyr)
library(vegan)
library(rcompanion)
library(ggplot2)
library(ggfortify)


#Load Drought Datasets
droughttraitsnoNA <- read.csv("Species_drought_traits_withoutNA.csv")
droughttraitsnoSLA <- read.csv("Species_droguht_traits_noSLA.csv")
droughttraitsnoRoot <- read.csv("Species_drought_traits_noroot.csv")

#Load Fire Datasets
firetraitsnoNA <- read.csv("Fixed_fire_traits_noNA.csv")
firetraitsnoLitter <- read.csv("Fixed_fire_traits_noLitter.csv")
firetraitsnoResprout <- read.csv("Fixed_fire_traits_noResprout.csv")
firetraitsnoSLA <- read.csv("Fixed_fire_traits_noSLA.csv")

#Check distributions and transform data if necessary
#Drought
plotNormalHistogram(droughttraitsnoNA$P50)
plotNormalHistogram(droughttraitsnoNA$Stem_conductivity)
plotNormalHistogram(droughttraitsnoNA$Root_depth)
plotNormalHistogram(droughttraitsnoNA$SLA)

#P50
hist(droughttraitsnoNA$P50) #No transform required

#Stem Conductivity
hist(droughttraitsnoNA$Stem_conductivity)
hist(log10(droughttraitsnoNA$Stem_conductivity)) #log10 transform seems to work the best
#Apply transform to all datsets
droughttraitsnoNA$Stem_conductivity <- log10(droughttraitsnoNA$Stem_conductivity)
droughttraitsnoRoot$Stem_conductivity <- log10(droughttraitsnoRoot$Stem_conductivity)
droughttraitsnoSLA$Stem_conductivity <- log10(droughttraitsnoSLA$Stem_conductivity)

#Root Depth
hist(droughttraitsnoNA$Root_depth)
hist(log(droughttraitsnoNA$Root_depth)) #log works the best
#Apply transform to all datasets
droughttraitsnoNA$Root_depth <- log(droughttraitsnoNA$Root_depth)
droughttraitsnoSLA$Root_depth <- log(droughttraitsnoSLA$Root_depth)

#SLA
hist(droughttraitsnoNA$SLA)
hist(log(droughttraitsnoNA$SLA)) #log works the best
#Apply transform to all datasets
droughttraitsnoNA$SLA <- log(droughttraitsnoNA$SLA)
droughttraitsnoRoot$SLA <- log(droughttraitsnoRoot$SLA)
#Remove NA's from No Root produced during log transformation
droughttraitsnoRoot <- na.omit(droughttraitsnoRoot)

#__________________________________________________________________________________________
#Fire
plotNormalHistogram(firetraitsnoNA$Bark_thickness)
plotNormalHistogram(firetraitsnoNA$SLA)
plotNormalHistogram(firetraitsnoNA$Litter_decomposition)
plotNormalHistogram(firetraitsnoNA$Resprouting)
plotNormalHistogram(firetraitsnoNA$Seed_mass)

#Bark Thickness
hist(firetraitsnoNA$Bark_thickness)
hist(log(firetraitsnoNA$Bark_thickness)) #log works the best
#Apply transforms to all datasets
firetraitsnoNA$Bark_thickness <- log(firetraitsnoNA$Bark_thickness)
firetraitsnoLitter$Bark_thickness <- log(firetraitsnoLitter$Bark_thickness)
firetraitsnoResprout$Bark_thickness <- log(firetraitsnoResprout$Bark_thickness)
firetraitsnoSLA$Bark_thickness <- log(firetraitsnoSLA$Bark_thickness)

#SLA
hist(firetraitsnoNA$SLA)
hist(log(firetraitsnoNA$SLA)) #log works the best
#Apply transforms to all datasets
firetraitsnoNA$SLA <- log(firetraitsnoNA$SLA)
firetraitsnoLitter$SLA <- log(firetraitsnoLitter$SLA)
firetraitsnoResprout$SLA <- log(firetraitsnoResprout$SLA)

# Litter Decomposition
hist(firetraitsnoNA$Litter_decomposition)
hist(sqrt(firetraitsnoNA$Litter_decomposition)) #sqrt works the best
#Apply transforms to all datasets
firetraitsnoNA$Litter_decomposition <- sqrt(firetraitsnoNA$Litter_decomposition)
firetraitsnoResprout$Litter_decomposition <- sqrt(firetraitsnoResprout$Litter_decomposition)
firetraitsnoSLA$Litter_decomposition <- sqrt(firetraitsnoSLA$Litter_decomposition)

#Resprouting Ability
hist(firetraitsnoNA$Resprouting)
#Categorical data so no point in considering the distribution or transforming

#Seed Dry Mass
hist(firetraitsnoNA$Seed_mass)
hist(log(firetraitsnoNA$Seed_mass)) #log works the best
#Apply transforms to all datasets
firetraitsnoNA$Seed_mass <- log(firetraitsnoNA$Seed_mass)
firetraitsnoLitter$Seed_mass <- log(firetraitsnoLitter$Seed_mass)
firetraitsnoResprout$Seed_mass <- log(firetraitsnoResprout$Seed_mass)
firetraitsnoSLA$Seed_mass <- log(firetraitsnoSLA$Seed_mass)

#Normality has been checked for all traits, and data have been transformed as necessary
#________________________________________________________________________________________
# Explore some relationships between traits

require(car)

# Litter decomposition rate and SLA for fire
litter <- firetraitsnoNA$Litter_decomposition
SLAf <- firetraitsnoNA$SLA

plot(litter, SLAf, main = "Litter Decomposition Rate vs SLA", sub = "Correlation Coefficient = 0.478",
     xlab = "Litter decomposition rate", ylab = "SLA", pch = 19, frame = TRUE)+
  abline(lm(SLAf~litter, data = firetraitsnoNA))

cor(litter, SLAf) #0.478 

# Stem conductivity and root depth for drought

stem <- droughttraitsnoNA$Stem_conductivity
root <- droughttraitsnoNA$Root_depth

plot(stem, root, main = "Stem Conductivity vs Rooting Depth", sub = "Correlation coefficient = -0.317",
     xlab = "Stem Conductivity", ylab = "Rooting Depth", pch = 19, frame = TRUE)+
  abline(lm(root~stem, data = droughttraitsnoNA))

cor(stem, root) # -0.317

#__________________________________________________________________________________________
#PCA: Drought Traits

#All Traits
all.drought.pca = prcomp(droughttraitsnoNA[2:5], center = T, scale = T)

sd.alldrought = all.drought.pca$sdev
loadings.alldrought = all.drought.pca$rotation
trait.names.alldrought = colnames(droughttraitsnoNA[2:5])
scores.alldrought = as.data.frame(all.drought.pca$x)
scores.alldrought$Species = droughttraitsnoNA$Species
summary(all.drought.pca)

# Save loadings for all drought traits
setwd("~/Dropbox/Trait_Data/Data/Trait_Data")
write.csv(loadings.alldrought, "PCA loadings all drought.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alldrought, "PCA scores all drought.csv")

#Biplot
autoplot(all.drought.pca, data = droughttraitsnoNA, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.5,0.6), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for All Drought Traits")+
  theme_classic()

#Save biplot
setwd("~/Dropbox/Trait_Data/Data/Trait_Data")
ggsave("all_drought_PCA.pdf", width = 6, height = 5, units = "in")

#----------------------------------------
#Drought no SLA
drought.noSLA.pca = prcomp(droughttraitsnoSLA[2:4], center = T, scale = T)

sd.noSLAd = drought.noSLA.pca$sdev
loadings.noSLAd = drought.noSLA.pca$rotation
trait.names.noSLAd = colnames(droughttraitsnoSLA[2:4])
scores.noSLAd = as.data.frame(drought.noSLA.pca$x)
scores.noSLAd$Species = droughttraitsnoSLA$Species
summary(drought.noSLA.pca)

# Save loadings for drought traits minus SLA
write.csv(loadings.noSLAd, "PCA loadings drought noSLA.csv", row.names = TRUE)

#Save species scores
write.csv(scores.noSLAd, "PCA scores drought noSLA.csv")

#Biplot
autoplot(drought.noSLA.pca, data = droughttraitsnoSLA, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.6), ylim = c(-0.7,0.5))+
  ggtitle(label = "PCA for Drought Traits Excluding SLA")+
  theme_classic()

#Save biplot
ggsave("drought_noSLA_PCA.pdf", width = 6, height = 5, units = "in")

#---------------------------------------
#Drought no Root Depth
drought.noRoot.pca = prcomp(droughttraitsnoRoot[2:4], center = T, scale = T)

sd.noRoot = drought.noRoot.pca$sdev
loadings.noRoot = drought.noRoot.pca$rotation
trait.names.noRoot = colnames(droughttraitsnoRoot[2:4])
scores.noRoot = as.data.frame(drought.noRoot.pca$x)
scores.noRoot$Species = droughttraitsnoRoot$Species
summary(drought.noRoot.pca)

# Save loadings for drought traits minus root depth
write.csv(loadings.noRoot, "PCA loadings drought noRoot.csv", row.names = TRUE)

#Save species scores
write.csv(scores.noRoot, "PCA scores drought noRoot.csv")

#Biplot
autoplot(drought.noRoot.pca, data = droughttraitsnoRoot, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3, geom = 'point', loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.6), ylim = c(-0.5,0.6))+
  ggtitle(label = "PCA for Drought Traits Excluding Root Depth")+
  theme_classic()

#Save biplot
ggsave("drought_noRoot_PCA.pdf", width = 6, height = 5, units = "in")

#____________________________________________________________________________________________
#PCA: Fire Traits

#All Traits
all.fire.pca = prcomp(firetraitsnoNA[2:6], center = T, scale = T)

sd.allfire = all.fire.pca$sdev
loadings.allfire = all.fire.pca$rotation
trait.names.allfire = colnames(firetraitsnoNA[2:6])
scores.allfire = as.data.frame(all.fire.pca$x)
scores.allfire$Species = firetraitsnoNA$Species
summary(all.fire.pca)

# Save loadings for all fire traits
write.csv(loadings.allfire, "PCA loadings all fire.csv", row.names = TRUE)

#Save species scores
write.csv(scores.allfire, "PCA scores all fire.csv")

#Biplot
autoplot(all.fire.pca, data = firetraitsnoNA, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3,loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.5,0.7), ylim = c(-0.5,0.7))+
  ggtitle(label = "PCA for All Fire Traits")+
  theme_classic()

#Save biplot
ggsave("all_fire_PCA.pdf", width = 6, height = 5, units = "in")

#----------------------------------------------
#Fire no Litter decomposition rate 
fire.nolitter.pca = prcomp(firetraitsnoLitter[2:5], center = T, scale = T)

sd.nolitter = fire.nolitter.pca$sdev
loadings.nolitter = fire.nolitter.pca$rotation
trait.names.nolitter = colnames(firetraitsnoLitter[2:5])
scores.nolitter = as.data.frame(fire.nolitter.pca$x)
scores.nolitter$Species = firetraitsnoLitter$Species
summary(fire.nolitter.pca)

# Save loadings for fire traits minus litter decomposition rate
write.csv(loadings.nolitter, "PCA loadings fire noLitter.csv", row.names = TRUE)

#Save species scores
write.csv(scores.nolitter, "PCA scores fire noLitter.csv")

#Biplot
autoplot(fire.nolitter.pca, data = firetraitsnoLitter, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3, loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.6,0.7), ylim = c(-0.5,0.7))+
  ggtitle(label = "PCA for Fire Traits Excluding Litter Decomposition Rate")+
  theme_classic()

#Save biplot
ggsave("fire_noLitter_PCA.pdf", width = 6, height = 5, units = "in")

#--------------------------------------------
#Fire no Resprout Ability
fire.noresprout.pca = prcomp(firetraitsnoResprout[2:5], center = T, scale = T)

sd.noresprout = fire.noresprout.pca$sdev
loadings.noresprout = fire.noresprout.pca$rotation
trait.names.noresprout = colnames(firetraitsnoResprout[2:5])
scores.noresprout = as.data.frame(fire.noresprout.pca$x)
scores.noresprout$Species = firetraitsnoResprout$Species
summary(fire.noresprout.pca)

# Save loadings for fire traits minus resprout ability
write.csv(loadings.noresprout, "PCA loadings fire noResprout.csv", row.names = TRUE)

#Save species scores
write.csv(scores.noresprout, "PCA scores fire noResprout.csv")

#Biplot
autoplot(fire.noresprout.pca, data = firetraitsnoResprout, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3, loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.6,0.7), ylim = c(-0.7,0.6))+
  ggtitle(label = "PCA for Fire Traits Excluding Resprout Ability")+
  theme_classic()

#Save biplot
ggsave("fire_noResprout_PCA.pdf", width = 6, height = 5, units = "in")

#------------------------------------------
#Fire no SLA
fire.noSLA.pca = prcomp(firetraitsnoSLA[2:5], center = T, scale = T)

sd.noSLAf = fire.noSLA.pca$sdev
loadings.noSLAf = fire.noSLA.pca$rotation
trait.names.noSLAf = colnames(firetraitsnoSLA[2:5])
scores.noSLAf = as.data.frame(fire.noSLA.pca$x)
scores.noSLAf$Species = firetraitsnoSLA$Species
summary(fire.noSLA.pca)

# Save loadings for fire traits minus SLA
write.csv(loadings.noSLAf, "PCA loadings fire noSLA.csv", row.names = TRUE)

#Save species scores
write.csv(scores.noSLAf, "PCA scores fire noSLA.csv")

#Biplot
autoplot(fire.noSLA.pca, data = firetraitsnoSLA, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3, loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.7,0.6), ylim = c(-0.6,0.5))+
  ggtitle(label = "PCA for Fire Traits Excluding SLA")+
  theme_classic()

#Save biplot
ggsave("fire_noSLA_PCA.pdf", width = 6, height = 5, units = "in")

#------------------------------------------
#Fire no Resprout and no Seed Mass

fire.noResproutSeed.pca = prcomp(firetraitsnoNA[2:4], center = T, scale = T)

sd.norespoutseed = fire.noResproutSeed.pca$sdev
loadings.noresproutseed = fire.noResproutSeed.pca$rotation
trait.names.norespoutseed = colnames(firetraitsnoNA[2:4])
scores.norespoutseed = as.data.frame(fire.noResproutSeed.pca$x)
scores.norespoutseed$Species = firetraitsnoNA$Species
summary(fire.noResproutSeed.pca)

# Save loadings for fire traits minus Resprout Ability and Seed Mass
write.csv(loadings.noresproutseed, "PCA loadings fire noResproutSeed.csv", row.names = TRUE)

#Save species scores
write.csv(scores.norespoutseed, "PCA scores fire noResproutSeed.csv")

#Biplot
autoplot(fire.noResproutSeed.pca, data = firetraitsnoNA, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3,loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.5,0.6), ylim = c(-0.8,0.6))+
  ggtitle(label = "PCA for Fire Traits Excluding Resprout Ability and Seed Mass")+
  theme_classic()

#Save biplot
ggsave("fire_noResproutSeed_PCA.pdf", width = 6, height = 5, units = "in")

#------------------------------------------
#Fire no Seed Mass

fire.noSeed.pca = prcomp(firetraitsnoNA[2:5], center = T, scale = T)

sd.noseed = fire.noSeed.pca$sdev
loadings.noseed = fire.noSeed.pca$rotation
trait.names.noseed = colnames(firetraitsnoNA[2:5])
scores.noseed = as.data.frame(fire.noSeed.pca$x)
scores.noseed$Species = firetraitsnoNA$Species
summary(fire.noSeed.pca)

# Save loadings for fire traits minus Resprout Ability and Seed Mass
write.csv(loadings.noseed, "PCA loadings fire noSeedMass.csv", row.names = TRUE)

#Save species scores
write.csv(scores.noseed, "PCA scores fire noSeedMass.csv")

#Biplot
autoplot(fire.noSeed.pca, data = firetraitsnoNA, label = T, check_overlap = T, label.label = "Species", 
         label.size = 3,loadings = T, loadings.colour = 'red',
         loadings.label = T, loadings.label.size = 3, loadings.label.repel = T, 
         xlim = c(-0.6,0.6), ylim = c(-0.5,0.7))+
  ggtitle(label = "PCA for Fire Traits Excluding Seed Mass")+
  theme_classic()

#Save biplot
ggsave("fire_noSeedMass_PCA.pdf", width = 6, height = 5, units = "in")


#--------------------------------------------------
#PCA of dataset with all traits (drought and fire)

#Used the SLA measurement from the drought dataset because that one was transformed,
#and the drought dataset had more species than the fire did
alltraits <- merge(droughttraitsnoNA, firetraitsnoSLA)

all.traits.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.alltraits = all.traits.pca$sdev
loadings.alltraits = all.traits.pca$rotation
trait.names.alltraits = colnames(alltraits[2:9])
scores.alltraits = as.data.frame(all.traits.pca$x)
scores.alltraits$Species = alltraits$Species
summary(all.traits.pca)

#Save loadings for all drought traits
setwd("~/Dropbox/Trait_Data/Data/Trait_Data")
write.csv(loadings.alltraits, "PCA loadings all traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "PCA scores all traits.csv")

#Biplot
autoplot(all.traits.pca, data = alltraits, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.9,0.5), ylim = c(-0.6,0.6))+
  ggtitle(label = "PCA for All Traits")+
  theme_classic()

#Save biplot
ggsave("all_traits_PCA.pdf", width = 6, height = 5, units = "in")

#----------------------------------------------
#Testing a couple other combinations 


#PCA of all traits excluding rooting depth and bark thickness

traitsnorootbark <- arrange(alltraits, Species)%>%
  select(Species, P50, Stem_conductivity, SLA, Litter_decomposition, Resprouting, Seed_mass)

traitsnorootbark.pca <- prcomp(traitsnorootbark[2:7], center = T, scale = T)

summary(traitsnorootbark.pca)

autoplot(traitsnorootbark.pca, data = traitsnorootbark, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.9,0.5), ylim = c(-0.6,0.6))+
  ggtitle(label = "PCA for Traits Excluding Rooting Depth and Bark Thickness")+
  theme_classic()

#PCA of all traits excluding rooting depth, bark thickness and resprout ability

traitsnoRBR <- arrange(alltraits, Species)%>%
  select(Species, P50, Stem_conductivity, SLA, Litter_decomposition, Seed_mass)

traitsnoRBR.pca <- prcomp(traitsnoRBR[2:6], center = T, scale = T)

summary(traitsnoRBR.pca)

autoplot(traitsnoRBR.pca, data = traitsnoRBR, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.7), ylim = c(-0.6,0.6))+
  ggtitle(label = "PCA for Traits Excluding Rooting Depth, Bark Thickness, and Resprout Ability")+
  theme_classic()

#PCA of all traits excluding seed mass, resprout ability 

traitsnoseedresprout <- arrange(alltraits, Species)%>%
  select(Species, P50, Root_depth, Stem_conductivity, SLA, Bark_thickness, Litter_decomposition)

traitsnoseedresprout.pca <- prcomp(traitsnoseedresprout[2:6], center = T, scale = T)

summary(traitsnoseedresprout.pca)

autoplot(traitsnoseedresprout.pca, data = traitsnoseedresprout, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.7), ylim = c(-0.6,0.6))+
  ggtitle(label = "PCA for Traits Excluding Seed Mass and Resprout Ability")+
  theme_classic()

#PCA of all traits excluding SLA and resprout ability

traitsnoSLAresprout <- arrange(alltraits, Species)%>%
  select(Species, P50, Root_depth, Stem_conductivity, Bark_thickness, Litter_decomposition, Seed_mass)

traitsnoSLAresprout.pca <- prcomp(traitsnoSLAresprout[2:6], center = T, scale = T)

summary(traitsnoSLAresprout.pca)

autoplot(traitsnoSLAresprout.pca, data = traitsnoSLAresprout, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.7), ylim = c(-0.6,0.6))+
  ggtitle(label = "PCA for Traits Excluding SLA and Resprout Ability")+
  theme_classic()

#PCA of all traits excluding litter decomp and seed mass 

traitsnolitterseed <- arrange(alltraits, Species)%>%
  select(Species, SLA, P50, Root_depth, Stem_conductivity, Bark_thickness, Resprouting)

traitsnolitterseed.pca <- prcomp(traitsnolitterseed[2:7], center = T, scale = T)

summary(traitsnolitterseed.pca)

autoplot(traitsnolitterseed.pca, data = traitsnolitterseed, label = T, check_overlap = T, label.label = "Species",
         label.size = 3, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.7,0.6), ylim = c(-0.8,0.5))+
  ggtitle(label = "PCA for Traits Excluding litter decomp and seed mass")+
  theme_classic()


#-----------------------------------------------------------------------
library(philentropy)

distance(droughttraitsnoNA[2:5], method = "euclidean")

