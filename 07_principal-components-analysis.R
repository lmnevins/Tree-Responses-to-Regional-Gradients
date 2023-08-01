######################################################################
#
#  Tree Responses to Regional Gradients -- PCA with imputed trait data
#
# L. McKinley Nevins, laura.nevins@wsu.edu, 22 July 2020
#
##      GNU General Public License, Version 3.0    ###################

# PCAs from the imputed dataset produced using the phylopars function in the Rphylopars package.

setwd("~/Dropbox/Trait_Data/Data/Nick_output")

library(dplyr)
require(tidyr)
library(vegan)
library(rcompanion)
require(ggplot2)
require(ggfortify)
library(oompaBase)
library(oompaData)
library(mclust)
library(ClassDiscovery)
library(kernlab)
library(changepoint)
library(cpm)
library(PCDimension)

# Load imputed data file of all traits 
alltraits <- read.csv("phy_imp_output.csv")


####################################
#PCA's
#All Traits
alltraits.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(alltraits[2:9])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$Species = alltraits$Species
summary(alltraits.pca)

# Save loadings for all traits
setwd("~/Dropbox/WSU/Cascade_Project/J of Ecol Submission/Resubmission2")
write.csv(loadings.alltraits, "PCA_loadings_alltraits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "PCA_scores_alltraits.csv")

#General Biplot to Assess Trends
plot <- ggplot2::autoplot(alltraits.pca, data = alltraits, label = T, check_overlap = T, label.label = "Species",
                          label.size = 3, loadings = T, loadings.colour = 'blue',
                          loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
                          xlim = c(-0.5,0.6), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for All Traits")+
  theme_classic()

print(plot)

##Broken-Stick test for the significance of the loadings
bs_pca <- alltraits.pca(t(ranData))
######################################################################
##Adding Species Codes to the dataset

Species_codes <- c('ABAM','ABCO','ABGR','ABLA','ABMA','ABPR','ACMA3','AECA','ALRH2','ALRU2','ARME','BEOC2','BEPA',
                   'CADE27','CHLA','CHNO','CHCHC4','CONU4','CUSA3','FRLA','JUOC','LALY','LAOC','LIDE3','MAFU','PIBR',
                   'PIEN','PISI','PIAL','PIAT','PICO','PIJE','PILA','PIMO3','PIMU','PIPO','PISA2','POTR5','POBAT','PREM',
                   'PRVI','PSME','QUAG','QUCH2','QUDO','QUGA4','QUKE','QULO','QUWI2','SESE3','TABR2','THPL','TOCA','TSHE','TSME','UMCA')

alltraits<- cbind(alltraits,Species_codes)

##Adding Hardwood or Conifer Designation to the dataset

Type <- c('C','C','C','C','C','C','H','H','H','H','H','H','H','C','C','C','H','H',
          'C','H','C','C','C','H','H','C','C','C','C','C','C','C','C','C','C','C',
          'C','H','H','H','H','C','H','H','H','H','H','H','H','C','C','C','C','C',
          'C','H')

alltraits<- cbind(alltraits,Type)


#PCA for All Traits Using Species Codes Instead

spcodes.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.spcodes = spcodes.pca$sdev
loadings.spcodes = spcodes.pca$rotation
trait.names.spcodes = colnames(alltraits[2:9])
scores.spcodes = as.data.frame(spcodes.pca$x)
scores.spcodes$Species_codes = alltraits$Species_codes
summary(spcodes.pca)

#Same as the loadings from above 
write.csv(loadings.spcodes, "PCA_loadings_spcodes.csv", row.names = TRUE)

#Set colors from conifers and hardwoods
pallette <- c("#00571A","#19ABE0")

#Biplot
plot <- autoplot(spcodes.pca, data = alltraits, colour = 'Type', label = T,
                 check_overlap = T, label.label = "Species_codes", repel = T,
                 label.size = 3, loadings = T, loadings.colour = 'black',
                 loadings.label = F, loadings.label.size = 0.0, loadings.label.colour = 'black', 
                 loadings.label.repel = T, xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  theme_classic()+
  geom_point() +
  scale_colour_manual(values=pallette, 
                      name="Tree Type",
                      breaks=c("C", "H"),
                      labels=c("Conifer", "Hardwood"))+
  theme(legend.position=c(0.12,0.16))+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))

plot

delete_layers(plot, "GeomPoint")

ggsave("PCA_coniferhardwood.pdf", width = 6, height = 5, units = "in")