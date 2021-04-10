# Creating plot for all traits from the phylopars data.

setwd("~/Dropbox/Trait_Data/Data/Nick_output")

library(dplyr)
require(tidyr)
library(vegan)
library(rcompanion)
require(ggplot2)
require(ggfortify)
library(ggrepel)
library(gginnards)

devtools::install_github("AliciaSchep/gglabeller") 

library(gglabeller)

# Load imputed data file of all traits 
alltraits <- read.csv("phy_imp_output.csv")

#PCA's
#All Traits
alltraits.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(alltraits[2:9])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$Species = alltraits$Species
summary(alltraits.pca)

#Biplot
autoplot(alltraits.pca, data = alltraits, label = T, check_overlap = T, label.label = "Species",
                  label.size = 2.5, loadings = T, loadings.colour = 'blue',
                  loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
                  xlim = c(-0.5,0.5), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for All Traits")+
  theme_classic()


#Save biplot
setwd("~/Dropbox/Trait_Data/Results/PCA results")
ggsave("PCA_phylopars_alltraits.pdf", width = 6, height = 5, units = "in")

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

##Adding abundance ranking to the dataset

Rank <- c('6','5','7','12','25','26','19','49','45','11','18','56','38','14','36','33','28','41','53','44','15','46',
          '16','8','52','55','17','29','32','43','4','21','24','27','50','2','42','37','39','40','51','1','48','13','30',
          '23','20','47','34','22','35','9','54','3','10','31')

alltraits<- cbind(alltraits,Rank)

alltraits<- mutate(alltraits$ranknames = dplyr::case_when(Rank > 0 ~ Species_codes, Rank < 21 ~ Species_codes, TRUE ~ ""))




#PCA for All Traits Using Species Codes Instead

spcodes.pca = prcomp(alltraits[2:9], center = T, scale = T)

sd.spcodes = spcodes.pca$sdev
loadings.spcodes = spcodes.pca$rotation
trait.names.spcodes = colnames(alltraits[2:9])
scores.spcodes = as.data.frame(spcodes.pca$x)
scores.spcodes$Species_codes = alltraits$Species_codes
summary(spcodes.pca)

#Set colors from conifers and hardwoods
pallette <- c("#00571A","#19ABE0")

#Biplot
 plot <- autoplot(spcodes.pca, data = alltraits, colour = 'Type', label = T,
         check_overlap = T, label.label = "Species_codes", repel = T,
         label.size = 3, loadings = T, loadings.colour = 'black',
         loadings.label = T, loadings.label.size = 3.5, loadings.label.colour = 'black', 
         loadings.label.repel = T, xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for All Traits using Species Codes")+
   geom_point() +
  theme_classic()+
  scale_colour_manual(values=pallette, 
                      name="Tree Type",
                      breaks=c("C", "H"),
                      labels=c("Conifer", "Hardwood"))+
  theme(legend.position=c(0.12,0.16))+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))
 
   delete_layers(plot, "GeomPoint")


ggsave("PCA_coniferhardwood.pdf", width = 6, height = 5, units = "in")



##Testing gglabeller function

plot <- autoplot(spcodes.pca, data = alltraits, colour = 'Type', label = F, 
                 loadings = T, loadings.colour = 'black',
                 loadings.label = T, loadings.label.size = 4, 
                 loadings.label.colour = 'black', 
                 loadings.label.repel = T, xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for All Traits using Species Codes")+
  theme_classic()+
  scale_colour_manual(values=pallette, 
                      name="Tree Type",
                      breaks=c("C", "H"),
                      labels=c("Conifer", "Hardwood"))+
  theme(legend.position=c(0.12,0.16))+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))


gglabeller_species <- gglabeller(plot, aes(label = Species_codes))




autoplot(spcodes.pca, data = alltraits, colour = 'Type', loadings = T, loadings.colour = 'black',
         loadings.label = T, loadings.label.size = 4, loadings.label.colour = 'black', 
         loadings.label.repel = T, xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for All Traits using Species Codes")+
  theme_classic()+
  scale_colour_manual(values=pallette, 
                      name="Tree Type",
                      breaks=c("C", "H"),
                      labels=c("Conifer", "Hardwood"))+
  theme(legend.position=c(0.12,0.16))+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))+
  geom_label_repel(data = dplyr::case_when(Rank > 0 ~ "Species_codes", Rank < 21 ~ "Species_codes", TRUE ~ ""))



##Get grayscale color scheme
#Set colors from conifers and hardwoods
pallette2 <- c("#999999","#000000")

#Biplot
autoplot(spcodes.pca, data = alltraits, colour = 'Type', label = T, 
         check_overlap = T, label.label = "Species_codes", repel = T,
         label.size = 2.5, loadings = T, loadings.colour = 'black',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'red', 
         loadings.label.repel = T, xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for All Traits using Species Codes")+
  theme_classic()+
  scale_colour_manual(values=pallette2, 
                      name="Tree Type",
                      breaks=c("C", "H"),
                      labels=c("Conifer", "Hardwood"))+
  theme(legend.position=c(0.12,0.16))+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="black", size = 12))





#----------------------------------------
####Other Explorations:
#All Traits Excluding Resprouting
traitsResprout.pca = prcomp(alltraits[2:8], center = T, scale = T)

sd.traitsResprout = traitsResprout.pca$sdev
loadings.traitsResprout = traitsResprout.pca$rotation
trait.names.traitsResprout = colnames(alltraits[2:8])
scores.traitsResprout = as.data.frame(traitsResprout.pca$x)
scores.traitsResprout$Species = alltraits$Species
summary(traitsResprout.pca)

#Biplot
autoplot(traitsResprout.pca, data = alltraits, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.4), ylim = c(-0.5,0.3))+
  ggtitle(label = "PCA for Traits Excluding Resprouting")+
  theme_classic()

#---------------------------------------------------------------------
#Separate angiosperms and gymnosperms 

angio <- slice(alltraits, 7:13,17,18,20,24,25,38:41,43:49,56 , preserve = FALSE)
gymno <- slice(alltraits, 1:6,14:16,19,21:23,26:37,42,50:55, preserve = FALSE)

#PCA for angiosperms
angio.pca = prcomp(angio[2:9], center = T, scale = T)

sd.angio = angio.pca$sdev
loadings.angio = angio.pca$rotation
trait.names.angio = colnames(angio[2:9])
scores.angio = as.data.frame(angio.pca$x)
scores.angio$Species = angio$Species
summary(angio.pca)

#Biplot
autoplot(angio.pca, data = angio, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.5,0.4), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for All Traits for Angiosperms")+
  theme_classic()


#Save biplot
ggsave("PCA_phylopars_angios.pdf", width = 6, height = 5, units = "in")

#------------------------------
#Just Drought for Angiosperms 
droughtangio <- select(angio,Species,SLA,P50,Stem_conductivity,Root_depth)

droughtangio.pca = prcomp(droughtangio[2:5], center = T, scale = T)

sd.droughtangio = droughtangio.pca$sdev
loadings.droughtangio = droughtangio.pca$rotation
trait.names.droughtangio = colnames(droughtangio[2:5])
scores.droughtangio = as.data.frame(droughtangio.pca$x)
scores.droughtangio$Species = droughtangio$Species
summary(droughtangio.pca)

#Biplot
autoplot(droughtangio.pca, data = droughtangio, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.6,0.5), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for Drought Traits for Angiosperms")+
  theme_classic()

#Just Fire for Angiosperms
fireangio <- select(angio,Species,Bark_thickness,Litter_decomposition,Seed_mass,Resprouting)

fireangio.pca = prcomp(fireangio[2:5], center = T, scale = T)

sd.fireangio = fireangio.pca$sdev
loadings.fireangio = fireangio.pca$rotation
trait.names.fireangio = colnames(fireangio[2:5])
scores.fireangio = as.data.frame(fireangio.pca$x)
scores.fireangio$Species = fireangio$Species
summary(fireangio.pca)

#Biplot
autoplot(fireangio.pca, data = fireangio, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.5), ylim = c(-0.7,0.5))+
  ggtitle(label = "PCA for Fire Traits for Angiosperms")+
  theme_classic()

#Fire for Angiosperms - No Resprouting
fireangioResprout <- select(angio,Species,Bark_thickness,Litter_decomposition,Seed_mass)

fireangioResprout.pca = prcomp(fireangioResprout[2:4], center = T, scale = T)

sd.fireangioResprout = fireangioResprout.pca$sdev
loadings.fireangioResprout = fireangioResprout.pca$rotation
trait.names.fireangioResprout = colnames(fireangioResprout[2:4])
scores.fireangioResprout = as.data.frame(fireangioResprout.pca$x)
scores.fireangioResprout$Species = fireangioResprout$Species
summary(fireangioResprout.pca)

#Biplot
autoplot(fireangioResprout.pca, data = fireangioResprout, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.5), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for Fire Traits for Angiosperms Excluding Resprouting")+
  theme_classic()


#=========================================================================
#PCA for Gymnosperms
gymno.pca = prcomp(gymno[2:9], center = T, scale = T)

sd.gymno = gymno.pca$sdev
loadings.gymno = gymno.pca$rotation
trait.names.gymno = colnames(gymno[2:9])
scores.gymno = as.data.frame(gymno.pca$x)
scores.gymno$Species = gymno$Species
summary(gymno.pca)

#Biplot
autoplot(gymno.pca, data = gymno, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.5,0.5), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for All Traits for Gymnosperms")+
  theme_classic()

#Save biplot
ggsave("PCA_phylopars_gymnos.pdf", width = 6, height = 5, units = "in")


#Remove resprouting as a consideration 
gymno2.pca = prcomp(gymno[2:8], center = T, scale = T)

sd.gymno2 = gymno2.pca$sdev
loadings.gymno2 = gymno2.pca$rotation
trait.names.gymno2 = colnames(gymno2[2:8])
scores.gymno2 = as.data.frame(gymno2.pca$x)
scores.gymno2$Species = gymno$Species
summary(gymno2.pca)

#Biplot
autoplot(gymno2.pca, data = gymno, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.5,0.5), ylim = c(-0.5,0.5))+
  ggtitle(label = "PCA for Gymnosperms Traits Excluding Resprouting")+
  theme_classic()

#------------------------------------
#Just Drought for Gymnosperms 
droughtgymno <- select(gymno,Species,SLA,P50,Stem_conductivity,Root_depth)

droughtgymno.pca = prcomp(droughtgymno[2:5], center = T, scale = T)

sd.droughtgymno = droughtgymno.pca$sdev
loadings.droughtgymno = droughtgymno.pca$rotation
trait.names.droughtgymno = colnames(droughtgymno[2:5])
scores.droughtgymno = as.data.frame(droughtgymno.pca$x)
scores.droughtgymno$Species = droughtgymno$Species
summary(droughtgymno.pca)

#Biplot
autoplot(droughtgymno.pca, data = droughtgymno, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.5), ylim = c(-0.4,0.5))+
  ggtitle(label = "PCA for Drought Traits for Gymnosperms")+
  theme_classic()

#Just Fire for Gymnosperms
firegymno <- select(gymno,Species,Bark_thickness,Litter_decomposition,Seed_mass,Resprouting)

firegymno.pca = prcomp(firegymno[2:5], center = T, scale = T)

sd.firegymno = firegymno.pca$sdev
loadings.firegymno = firegymno.pca$rotation
trait.names.firegymno = colnames(firegymno[2:5])
scores.firegymno = as.data.frame(firegymno.pca$x)
scores.firegymno$Species = firegymno$Species
summary(firegymno.pca)

#Biplot
autoplot(firegymno.pca, data = firegymno, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.5), ylim = c(-0.6,0.4))+
  ggtitle(label = "PCA for Fire Traits for Gymnosperms")+
  theme_classic()

#Fire for Gymnosperms - No Resprouting
firegymnoResprout <- select(gymno,Species,Bark_thickness,Litter_decomposition,Seed_mass)

firegymnoResprout.pca = prcomp(firegymnoResprout[2:4], center = T, scale = T)

sd.firegymnoResprout = firegymnoResprout.pca$sdev
loadings.firegymnoResprout = firegymnoResprout.pca$rotation
trait.names.firegymnoResprout = colnames(firegymnoResprout[2:4])
scores.firegymnoResprout = as.data.frame(firegymnoResprout.pca$x)
scores.firegymnoResprout$Species = firegymnoResprout$Species
summary(firegymnoResprout.pca)

#Biplot
autoplot(firegymnoResprout.pca, data = firegymnoResprout, label = T, check_overlap = T, label.label = "Species",
         label.size = 2.5, loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 3, loadings.label.colour = 'blue', loadings.label.repel = T, 
         xlim = c(-0.4,0.4), ylim = c(-0.4,0.7))+
  ggtitle(label = "PCA for Fire Traits for Gymnosperms Excluding Resprouting")+
  theme_classic()
