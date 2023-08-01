#############################################################################################################
#
#   Tree Responses to Regional Gradients -- Bootstrapping growth mixed-effects models & analysis of results 
#
#   L. McKinley Nevins, laura.nevins@wsu.edu, 6 June 2023
#   Jenny Zambrano, jenny.zambrano@wsu.edu 
#   Rob Smith, phytomosaic@gmail.com
#   
##      GNU General Public License, Version 3.0    #########################################################

require(ecole)     
require(data.table)
require(dplyr)
require(ggplot2)
require(lme4)
require(tidyverse)
require(car)
require(ggpubr)

####################################################################################

##load data
setwd("~/Documents/McKinley/RStudio/Data")
growth <- read.csv("growth_full.csv")

#################################################################################################################
##Generalized Linear Mixed Models 

#perform once to get coefficients 

##fire####
start_time <- Sys.time()
m_rgr_fire <- lme4::lmer(rgr ~ firelog * pc1              # second-order interaction term
                         + (1 + firelog * pc1 | spp)      # separate slopes for species 
                         + (dia_y1 + I(dia_y1^2))                  # standardized tree size 
                         + (1 | census_int)                          # interval between censuses 
                         + (1 | censusyr)                       # separate means for year
                         + (1 | plt_cn),                        # separate means for plot
                         data      = growth,                  
                         na.action = 'na.fail',                   # NA handling
                         control=lmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=2e6)))    

end_time <- Sys.time()
end_time - start_time # time elapsed ~ 3.8 minutes for growth unscaled on lab iMac

#check for singularity
isSingular(m_rgr_fire, tol = 1e-4)


#summaries
summary(m_rgr_fire)
coef(m_rgr_fire)
car::Anova(m_rgr_fire, type=3)
performance::r2_nakagawa(m_rgr_fire,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)


##cmd#####
start_time <- Sys.time()
m_rgr_cmd <- lme4::lmer(rgr ~ cmdlog * pc1              # second-order interaction term
                        + (1 + cmdlog * pc1 | spp)      # separate slopes for species 
                        + (dia_y1 + I(dia_y1^2))               # standardized tree size                            
                        + (1 | census_int)                          # interval between censuses 
                        + (1 | censusyr)                       # separate means for year
                        + (1 | plt_cn),                        # separate means for plot
                        data      = growth,                  
                        na.action = 'na.fail',                   # NA handling
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e6)))   
end_time <- Sys.time()
end_time - start_time # time elapsed ~ 2.77 minutes for growth unscaled 

#check for singularity
isSingular(m_rgr_cmd, tol = 1e-4)

#summaries
summary(m_rgr_cmd)
coef(m_rgr_cmd)
car::Anova(m_rgr_cmd, type=3)
performance::r2_nakagawa(m_rgr_cmd,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)



######### Bootstrapping ############################################################

##Growth - fire x PC1##################

#takes ~4 days to run 
nboots <- 999
cxs <- coef(m_rgr_fire)$spp[,c('firelog','pc1','firelog:pc1', 'I(dia_y1^2)')] # per-species coefficients for two main effects and interaction
a <- array(NA, dim=c(dim(cxs),nboots), dimnames=list(rownames(cxs),colnames(cxs),NULL))
# ! ! ! ! TIMEWARN ! ! ! ! about 3 min per iteration ! ! ! ! ! ! ! ! ! ! ! ! ! !
for (b in 1:nboots) {
  cat('iteration', b, 'of', nboots, '....... \n')
  nr <- NROW(growth)
  dx <- growth[sample(1:nr, size=nr, replace=FALSE),,] # draw with replacement
  m_rgr_fire <- lme4::lmer(rgr ~ firelog * pc1              # second-order interaction term
                           + (1 + firelog * pc1 | spp)      # separate slopes for species 
                           + (dia_y1 + I(dia_y1^2))               # standardized tree size 
                           + census_int                           # interval between censuses 
                           + (1 | censusyr)                       # separate means for year
                           + (1 | plt_cn),                        # separate means for plot
                           data      = dx,                  
                           na.action = 'na.fail',                   # NA handling
                           control=lmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e6)))  
  a[,,b] <- as.matrix(coef(m_rgr_fire)$spp[,c('firelog','pc1','firelog:pc1', 'I(dia_y1^2)')])
}

### move results array into clearer named file  
growth_fire_mod999 <- a                                 # resulting coefficients array in full

setwd("~/Documents/McKinley/Cascade_Project")
save(growth_fire_mod999, file="growth_fire_mod999.rda")



###################Growth - CMD x PC1 ##################

#takes ~3 days to run 
nboots <- 999
cxs <- coef(m_rgr_cmd)$spp[,c('cmdlog','pc1','cmdlog:pc1' ,'I(dia_y1^2)')] # per-species coefficients for two main effects, interaction, and diameter
a <- array(NA, dim=c(dim(cxs),nboots), dimnames=list(rownames(cxs),colnames(cxs),NULL))
# ! ! ! ! TIMEWARN ! ! ! ! about 3 min per iteration ! ! ! ! ! ! ! ! ! ! ! ! ! !
for (b in 1:nboots) {
  cat('iteration', b, 'of', nboots, '....... \n')
  nr <- NROW(growth)
  dx <- growth[sample(1:nr, size=nr, replace=FALSE),,] # draw with replacement
  start_time <- Sys.time()
  m_rgr_cmd <- lme4::lmer(rgr ~ cmdlog * pc1              # second-order interaction term
                          + (1 + cmdlog * pc1 | spp)      # separate slopes for species 
                          + (dia_y1 + I(dia_y1^2))               # standardized tree size                        
                          + census_int                           # interval between censuses 
                          + (1 | censusyr)                       # separate means for year
                          + (1 | plt_cn),                        # separate means for plot
                          data      = dx,                  
                          na.action = 'na.fail',                   # NA handling
                          control=lmerControl(optimizer="bobyqa",
                                              optCtrl=list(maxfun=2e5)))   
  end_time <- Sys.time()
  end_time - start_time
  a[,,b] <- as.matrix(coef(m_rgr_cmd)$spp[,c('cmdlog','pc1','cmdlog:pc1','I(dia_y1^2)')])
}


### move results array into clearer named file  
growth_cmd_mod999 <- a                                 # resulting coefficients array in full

setwd("~/Documents/McKinley/Cascade_Project")
save(growth_cmd_mod999, file="growth_cmd_mod999.rda")


####################################################################################################
##Load bootstrapped results arrays 

setwd("~/Documents/McKinley/Cascade_Project/")

#fire model 
load(file = 'growth_fire_mod999.rda')

#cmd model 
load(file = 'growth_cmd_mod999.rda')

#########Generating Final Plots from the Results####################################

####### Fire Results###########

### reshape to get mean coefficients
rgr_fire_allspp <- growth_fire_mod999          # resulting coefficients array in full
# reshape it
rgr_fire_allspp <- matrix(aperm(rgr_fire_allspp, c(1, 3, 2)), nrow = dim(rgr_fire_allspp)[1] * dim(rgr_fire_allspp)[3], dimnames = list(rep(rownames(rgr_fire_allspp),dim(rgr_fire_allspp)[3]), colnames(rgr_fire_allspp)))
rgr_fire_allspp <- data.frame(spp = rownames(rgr_fire_allspp), rgr_fire_allspp, row.names=NULL)
rgr_fire_allspp <- rgr_fire_allspp[order(rgr_fire_allspp$spp),,]

#54,945 lines in the datafile, because there were 55 species iterated 999 times 

#get mean, sd, and se for values across all species
rgr_fire_allspp_vals <- rgr_fire_allspp %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, se = ~sd(.x)/sqrt(length(.x)))))

####Calculate confidence intervals  

#for fire
firehist <- hist(rgr_fire_allspp$firelog)
quantile(rgr_fire_allspp$firelog, c(0.025, 0.975))

#set lower and upper bound from output 
firelog_low_CI <- -0.0018910661
firelog_upper_CI <- 0.0005967524

#as dataframes
firelog_low_CI <- as.data.frame(firelog_low_CI)
firelog_upper_CI <- as.data.frame(firelog_upper_CI)


##for pc1
pc1hist <- hist(rgr_fire_allspp$pc1)
quantile(rgr_fire_allspp$pc1, c(0.025, 0.975))

# set lower bound and upper bound
pc1_low_CI <- -0.000934641
pc1_upper_CI <- 0.001000648

#as dataframes
pc1_low_CI <- as.data.frame(pc1_low_CI)
pc1_upper_CI <- as.data.frame(pc1_upper_CI)


##for firelog x pc1
fire.pc1hist <- hist(rgr_fire_allspp$firelog.pc1)
quantile(rgr_fire_allspp$firelog.pc1, c(0.025, 0.975))

#set lower bound and upper bound
firelog.pc1_low_CI <- -0.0007518506
firelog.pc1_upper_CI <- 0.0006467613

#as dataframes
firelog.pc1_low_CI <- as.data.frame(firelog.pc1_low_CI)
firelog.pc1_upper_CI <- as.data.frame(firelog.pc1_upper_CI)

##for diameter
diameterhist <- hist(rgr_fire_allspp$I.dia_y1.2.)
quantile(rgr_fire_allspp$I.dia_y1.2., c(0.025, 0.975))

# Calculating lower bound and upper bound
dia_low_CI <- 0.0008109076
dia_upper_CI <- 0.0008109076

#as dataframes
dia_low_CI <- as.data.frame(dia_low_CI)
dia_upper_CI <- as.data.frame(dia_upper_CI)


#subset rgr_fire_allspp_vals to get the means 
rgr_fire_allspp_vals <- rgr_fire_allspp_vals %>%
  dplyr::select(firelog_mean, pc1_mean, firelog.pc1_mean, I.dia_y1.2._mean)


#bind everything together 
allspp_rgr_fire_intervals <- bind_cols(rgr_fire_allspp_vals, firelog_low_CI, firelog_upper_CI, 
                                       pc1_low_CI, pc1_upper_CI,
                                       firelog.pc1_low_CI, firelog.pc1_upper_CI,
                                       dia_low_CI, dia_upper_CI)

#combine into a new dataframe 
rgr_fire_allspp_vals <- data.frame(coefficient = c('Fire', 'Funct Neigh', 'Fire x Funct Neigh', 'Initial Tree Size'),
                                   mean = c(-0.0004050744, 9.104824e-05, -6.564223e-05, 0.0008109076),
                                   low_confint = c(-0.001891066, -0.000934641, -0.0007518506, 0.0008109076),
                                   upper_confint = c(0.0005967524, 0.001000648, 0.0006467613, 0.0008109076),
                                   sig = c('no', 'no', 'no', 'yes'))

rgr_fire_allspp_vals$sig <- as.factor(rgr_fire_allspp_vals$sig)

#plot coefficients and confidence intervals for all species combined 
allplot_rgrfire <-ggplot(rgr_fire_allspp_vals,aes(x=mean,y=coefficient)) +
  geom_point() + 
  geom_errorbar(aes(xmin=low_confint, xmax=upper_confint), colour="black", width=.1) +
  geom_point(color='black', shape=21, size=2, aes(fill=factor(sig)), show.legend = FALSE) + 
  scale_fill_manual(values=c('white', 'black')) +
  geom_line() +
  geom_vline(xintercept=0, size=1) +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.002, 0.00115) + #put on same scale as the CMD plot 
  ggtitle("") +
  theme_bw()


allplot_rgrfire + theme(axis.text.y = element_text(size = 10), 
                        axis.text.x = element_text(size = 10),)


######## assess each species individually #####################################################################################

### reshape the resulting bootstrapped coefficients array
rgr_fire_allspp <- growth_fire_mod999                  # resulting coefficients array in full
# reshape it
rgr_fire_15spp <- matrix(aperm(rgr_fire_allspp, c(1, 3, 2)), nrow = dim(rgr_fire_allspp)[1] * dim(rgr_fire_allspp)[3], dimnames = list(rep(rownames(rgr_fire_allspp),dim(rgr_fire_allspp)[3]), colnames(rgr_fire_allspp)))
rgr_fire_15spp <- data.frame(spp = rownames(rgr_fire_15spp), rgr_fire_15spp, row.names=NULL)
rgr_fire_15spp <- rgr_fire_15spp[order(rgr_fire_15spp$spp),,]

#subset 15 most abundant species 
rgr_fire_15spp = rgr_fire_15spp[rgr_fire_15spp$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                                          "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                                          "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 
#reduces to 14,985 observations for the 5 variables 

#get mean for each species 
rgr_fire_stats <- rgr_fire_15spp %>% 
  group_by(spp) %>% 
  summarise(across(everything(), list(mean = mean)))

#check the distributions of the coefficient values 
hist(rgr_fire_15spp$firelog)
hist(rgr_fire_15spp$pc1)
hist(rgr_fire_15spp$firelog.pc1)

#this calculates the upper and lower confidence intervals (here as the 2.5% and 97.5% quantiles, and organizes them by species)

#for fire
firelog_quants <- rgr_fire_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(firelog, c(0.025, 0.975)), "quantile", "firelog"))

#for pc1
pc1_quants <- rgr_fire_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(pc1, c(0.025, 0.975)), "quantile", "pc1"))

#for fire x pc1
firelog.pc1_quants <- rgr_fire_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(firelog.pc1, c(0.025, 0.975)), "quantile", "firelog.pc1"))


#make quantile a factor with two levels in each of the dataframes 
firelog_quants$quantile <- as.factor(firelog_quants$quantile)
pc1_quants$quantile <- as.factor(pc1_quants$quantile)
firelog.pc1_quants$quantile <- as.factor(firelog.pc1_quants$quantile)

#pull lower quantile out of each dataframe 
firelog_lowquant = firelog_quants[firelog_quants$quantile == "2.5%", ]
pc1_lowquant = pc1_quants[pc1_quants$quantile == "2.5%", ]
firelog.pc1_lowquant = firelog.pc1_quants[firelog.pc1_quants$quantile == "2.5%", ]

#bind together and rename some columns 
low_quants <- bind_cols(firelog_lowquant, pc1_lowquant, firelog.pc1_lowquant)
low_quants <- low_quants %>% dplyr::select(spp...1, firelog, pc1, firelog.pc1) %>% 
              rename(spp = spp...1, Fire = firelog, Funct_Neigh = pc1, FirexFunct_Neigh = firelog.pc1)

#pull upper quantile out of each dataframe 
firelog_upperquant = firelog_quants[firelog_quants$quantile == "97.5%", ]
pc1_upperquant = pc1_quants[pc1_quants$quantile == "97.5%", ]
firelog.pc1_upperquant = firelog.pc1_quants[firelog.pc1_quants$quantile == "97.5%", ]

#bind together and rename some columns 
upper_quants <- bind_cols(firelog_upperquant, pc1_upperquant, firelog.pc1_upperquant)
upper_quants <- upper_quants %>% dplyr::select(spp...1, firelog, pc1, firelog.pc1) %>% 
  rename(spp = spp...1, Fire = firelog, Funct_Neigh = pc1, FirexFunct_Neigh = firelog.pc1)


##bind together dataframes for means and lower and upper confidence intervals 

#subset means to reshape 
rgr_fire_15spp_means <- dplyr::select(.data = rgr_fire_stats, spp, firelog_mean, pc1_mean, firelog.pc1_mean)

rgr_fire_15spp_means <- rgr_fire_15spp_means %>% 
  rename(Fire = firelog_mean, Funct_Neigh = pc1_mean, FirexFunct_Neigh = firelog.pc1_mean)

#reshape the means dataset 
rgr_fire_15spp_means <- data.frame(spp = rownames(rgr_fire_15spp_means), rgr_fire_15spp_means, row.names=NULL)
rgr_fire_15spp_means <- rgr_fire_15spp_means[order(rgr_fire_15spp_means$spp),,]
rgr_fire_15spp_means  <- reshape(rgr_fire_15spp_means, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='mean',
                                 times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(rgr_fire_15spp_means)))

#reshape the low CI dataset 
rgr_fire_15spp_lowCI <- data.frame(spp = rownames(low_quants), low_quants, row.names=NULL)
rgr_fire_15spp_lowCI <- rgr_fire_15spp_lowCI[order(rgr_fire_15spp_lowCI$spp),,]
rgr_fire_15spp_lowCI  <- reshape(rgr_fire_15spp_lowCI, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='low_CI',
                                 times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(rgr_fire_15spp_lowCI)))

#reshape the upper CI dataset 
rgr_fire_15spp_upperCI <- data.frame(spp = rownames(upper_quants), upper_quants, row.names=NULL)
rgr_fire_15spp_upperCI <- rgr_fire_15spp_upperCI[order(rgr_fire_15spp_upperCI$spp),,]
rgr_fire_15spp_upperCI  <- reshape(rgr_fire_15spp_upperCI, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='upper_CI',
                                   times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                   new.row.names=1:prod(dim(rgr_fire_15spp_upperCI)))

#merge dataframes to get means and CI's all together 
rgr_fire_15spp_allses <- bind_cols(rgr_fire_15spp_means, rgr_fire_15spp_lowCI, rgr_fire_15spp_upperCI)
#select subset of columns and rename because the naming changed during merging 
#there's probably a cleaner way to do this but I don't know it
rgr_fire_15spp_allses <- dplyr::select(rgr_fire_15spp_allses, spp...1, spp.1...2, predictor...3, id...5, mean, low_CI, upper_CI) %>% 
  rename(rep = spp...1, spp = spp.1...2, id = id...5, predictor = predictor...3)

#reorder the species according to relative frequency 
#this doesn't objectively look like it changes anything, but it does plot in the right order
rgr_fire_15spp_allses$spp <- factor(rgr_fire_15spp_allses$spp, levels=rev(c("Pseudotsuga menziesii", "Pinus ponderosa",
                                                                            "Tsuga heterophylla", "Pinus contorta", "Abies concolor", "Abies amabilis",
                                                                            "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana",
                                                                            "Thuja plicata", "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis",
                                                                            "Calocedrus decurrens", "Juniperus occidentalis"))) # plot species by rel freq order

### plot the bootstrapped coefficients per species
rgr_fire_15 <- ggplot(data = rgr_fire_15spp_allses, aes(y = spp, x = mean, xmin=low_CI, xmax=upper_CI)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.5, color= "black", lwd=0.5) +
  geom_vline(xintercept=0, size=0.5, linetype='dashed') +
  facet_grid( ~ predictor, scales = 'fixed') +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-9e-04, 6e-04) +
  ggtitle("") +
  theme_bw()

rgr_fire_15 + theme(axis.text.y = element_text(face="italic", size = 12)) + theme(axis.text.x = element_text(size = 12)) + theme(strip.text.x = element_text(size = 12))
 

###########################################################################################################################

####### CMD Results###########

### reshape to get mean coefficients for each species 
rgr_cmd_allspp <- growth_cmd_mod999          # resulting coefficients array in full
# reshape it
rgr_cmd_allspp <- matrix(aperm(rgr_cmd_allspp, c(1, 3, 2)), nrow = dim(rgr_cmd_allspp)[1] * dim(rgr_cmd_allspp)[3], dimnames = list(rep(rownames(rgr_cmd_allspp),dim(rgr_cmd_allspp)[3]), colnames(rgr_cmd_allspp)))
rgr_cmd_allspp <- data.frame(spp = rownames(rgr_cmd_allspp), rgr_cmd_allspp, row.names=NULL)
rgr_cmd_allspp <- rgr_cmd_allspp[order(rgr_cmd_allspp$spp),,]

#54,945 lines in the datafile, because there were 55 species iterated 999 times 

#get mean, sd, and se for values across all species
rgr_cmd_allspp_vals <- rgr_cmd_allspp %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, se = ~sd(.x)/sqrt(length(.x)))))

####Calculating confidence intervals

#for cmd
cmdhist <- hist(rgr_cmd_allspp$cmdlog)
quantile(rgr_cmd_allspp$cmdlog, c(0.025, 0.975))

#set lower and upper bound from output 
cmdlog_low_CI <- -0.0009553685
cmdlog_upper_CI <- 0.0004624294

#as dataframes
cmdlog_low_CI <- as.data.frame(cmdlog_low_CI)
cmdlog_upper_CI <- as.data.frame(cmdlog_upper_CI)


##for pc1
pc1.2hist <- hist(rgr_cmd_allspp$pc1)
quantile(rgr_cmd_allspp$pc1, c(0.025, 0.975))

# set lower bound and upper bound
pc1_low_CI <- -0.0009755448
pc1_upper_CI <- 0.0006517349

#as dataframes
pc1_low_CI <- as.data.frame(pc1_low_CI)
pc1_upper_CI <- as.data.frame(pc1_upper_CI)


##for cmdlog x pc1
cmd.pc1hist <- hist(rgr_cmd_allspp$cmdlog.pc1)
quantile(rgr_cmd_allspp$cmdlog.pc1, c(0.025, 0.975))

#set lower bound and upper bound
cmdlog.pc1_low_CI <- -0.0006153433
cmdlog.pc1_upper_CI <- 0.0007477628

#as dataframes
cmdlog.pc1_low_CI <- as.data.frame(cmdlog.pc1_low_CI)
cmdlog.pc1_upper_CI <- as.data.frame(cmdlog.pc1_upper_CI)


##for diameter
diameterhist <- hist(rgr_cmd_allspp$I.dia_y1.2.)
quantile(rgr_cmd_allspp$I.dia_y1.2., c(0.025, 0.975))

# Calculating lower bound and upper bound
dia_low_CI <- 0.0008111933
dia_upper_CI <- 0.0008111933

#as dataframes
dia_low_CI <- as.data.frame(dia_low_CI)
dia_upper_CI <- as.data.frame(dia_upper_CI)


#subset rgr_cmd_allspp_vals to get the means 
rgr_cmd_allspp_vals <- rgr_cmd_allspp_vals %>%
  dplyr::select(cmdlog_mean, pc1_mean, cmdlog.pc1_mean, I.dia_y1.2._mean)


#bind everything together 
allspp_rgr_cmd_intervals <- bind_cols(rgr_cmd_allspp_vals, cmdlog_low_CI, cmdlog_upper_CI, 
                                       pc1_low_CI, pc1_upper_CI,
                                       cmdlog.pc1_low_CI, cmdlog.pc1_upper_CI,
                                       dia_low_CI, dia_upper_CI)

#combine into a new dataframe 
rgr_cmd_allspp_vals <- data.frame(coefficient = c('CMD', 'Funct Neigh', 'CMD x Funct Neigh', 'Initial Tree Size'),
                                   mean = c(-7.030465e-05, -8.270143e-05, 7.02099e-05, 0.0008111933),
                                   low_confint = c(-0.0009553685, -0.0009755448, -0.0006153433, 0.0008111933),
                                   upper_confint = c(0.0004624294, 0.0006517349, 0.0007477628, 0.0008111933),
                                  sig = c('no', 'no', 'no', 'yes'))


#plot coefficients and confidence intervals for all species combined 
allplot_rgrcmd <-ggplot(rgr_cmd_allspp_vals,aes(x=mean,y=coefficient)) +
  geom_point() + 
  geom_errorbar(aes(xmin=low_confint, xmax=upper_confint), colour="black", width=.1) +
  geom_point(color='black', shape=21, size=2, aes(fill=factor(sig)), show.legend = FALSE) + 
  scale_fill_manual(values=c('white', 'black')) +
  geom_line() +
  geom_vline(xintercept=0, size=1) +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.002, 0.00115) + #put on same scale as the fire plot 
  ggtitle("") +
  theme_bw()


allplot_rgrcmd + theme(axis.text.y = element_text(size = 10), 
                        axis.text.x = element_text(size = 10),)

######## assess each species individually #####################################################################################

### reshape the resulting bootstrapped coefficients array
rgr_cmd_allspp <- growth_cmd_mod999                  # resulting coefficients array in full
# reshape it
rgr_cmd_15spp <- matrix(aperm(rgr_cmd_allspp, c(1, 3, 2)), nrow = dim(rgr_cmd_allspp)[1] * dim(rgr_cmd_allspp)[3], dimnames = list(rep(rownames(rgr_cmd_allspp),dim(rgr_cmd_allspp)[3]), colnames(rgr_cmd_allspp)))
rgr_cmd_15spp <- data.frame(spp = rownames(rgr_cmd_15spp), rgr_cmd_15spp, row.names=NULL)
rgr_cmd_15spp <- rgr_cmd_15spp[order(rgr_cmd_15spp$spp),,]

#subset 15 most abundant species 
rgr_cmd_15spp = rgr_cmd_15spp[rgr_cmd_15spp$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                                          "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                                          "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 
#reduces to 14,985 observations for the 5 variables 

#get mean for each species 
rgr_cmd_stats <- rgr_cmd_15spp %>% 
  group_by(spp) %>% 
  summarise(across(everything(), list(mean = mean)))

#check the distributions of the coefficient values 
hist(rgr_cmd_15spp$cmdlog)
hist(rgr_cmd_15spp$pc1)
hist(rgr_cmd_15spp$cmdlog.pc1)

#this calculates the upper and lower confidence intervals (here as the 2.5% and 97.5% quantiles, and organizes them by species)

#for cmd
cmdlog_quants <- rgr_cmd_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(cmdlog, c(0.025, 0.975)), "quantile", "cmdlog"))

#for pc1
pc1_quants <- rgr_cmd_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(pc1, c(0.025, 0.975)), "quantile", "pc1"))

#for cmd x pc1
cmdlog.pc1_quants <- rgr_cmd_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(cmdlog.pc1, c(0.025, 0.975)), "quantile", "cmdlog.pc1"))


#make quantile a factor with two levels in each of the dataframes 
cmdlog_quants$quantile <- as.factor(cmdlog_quants$quantile)
pc1_quants$quantile <- as.factor(pc1_quants$quantile)
cmdlog.pc1_quants$quantile <- as.factor(cmdlog.pc1_quants$quantile)

#pull lower quantile out of each dataframe 
cmdlog_lowquant = cmdlog_quants[cmdlog_quants$quantile == "2.5%", ]
pc1_lowquant = pc1_quants[pc1_quants$quantile == "2.5%", ]
cmdlog.pc1_lowquant = cmdlog.pc1_quants[cmdlog.pc1_quants$quantile == "2.5%", ]

#bind together and rename some columns 
low_quants <- bind_cols(cmdlog_lowquant, pc1_lowquant, cmdlog.pc1_lowquant)
low_quants <- low_quants %>% dplyr::select(spp...1, cmdlog, pc1, cmdlog.pc1) %>% 
  rename(spp = spp...1, CMD = cmdlog, Funct_Neigh = pc1, CMDxFunct_Neigh = cmdlog.pc1)

#pull upper quantile out of each dataframe 
cmdlog_upperquant = cmdlog_quants[cmdlog_quants$quantile == "97.5%", ]
pc1_upperquant = pc1_quants[pc1_quants$quantile == "97.5%", ]
cmdlog.pc1_upperquant = cmdlog.pc1_quants[cmdlog.pc1_quants$quantile == "97.5%", ]

#bind together and rename some columns 
upper_quants <- bind_cols(cmdlog_upperquant, pc1_upperquant, cmdlog.pc1_upperquant)
upper_quants <- upper_quants %>% dplyr::select(spp...1, cmdlog, pc1, cmdlog.pc1) %>% 
  rename(spp = spp...1, CMD = cmdlog, Funct_Neigh = pc1, CMDxFunct_Neigh = cmdlog.pc1)


##bind together dataframes for means and lower and upper confidence intervals 
#subset means to reshape 
rgr_cmd_15spp_means <- dplyr::select(.data = rgr_cmd_stats, spp, cmdlog_mean, pc1_mean, cmdlog.pc1_mean)

rgr_cmd_15spp_means <- rgr_cmd_15spp_means %>% 
  rename(CMD = cmdlog_mean, Funct_Neigh = pc1_mean, CMDxFunct_Neigh = cmdlog.pc1_mean)

#reshape the means dataset 
rgr_cmd_15spp_means <- data.frame(spp = rownames(rgr_cmd_15spp_means), rgr_cmd_15spp_means, row.names=NULL)
rgr_cmd_15spp_means <- rgr_cmd_15spp_means[order(rgr_cmd_15spp_means$spp),,]
rgr_cmd_15spp_means  <- reshape(rgr_cmd_15spp_means, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='mean',
                                 times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(rgr_cmd_15spp_means)))

#reshape the low CI dataset 
rgr_cmd_15spp_lowCI <- data.frame(spp = rownames(low_quants), low_quants, row.names=NULL)
rgr_cmd_15spp_lowCI <- rgr_cmd_15spp_lowCI[order(rgr_cmd_15spp_lowCI$spp),,]
rgr_cmd_15spp_lowCI  <- reshape(rgr_cmd_15spp_lowCI, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='low_CI',
                                 times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(rgr_cmd_15spp_lowCI)))

#reshape the upper CI dataset 
rgr_cmd_15spp_upperCI <- data.frame(spp = rownames(upper_quants), upper_quants, row.names=NULL)
rgr_cmd_15spp_upperCI <- rgr_cmd_15spp_upperCI[order(rgr_cmd_15spp_upperCI$spp),,]
rgr_cmd_15spp_upperCI  <- reshape(rgr_cmd_15spp_upperCI, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='upper_CI',
                                   times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                   new.row.names=1:prod(dim(rgr_cmd_15spp_upperCI)))


#merge dataframes to get means and CI's all together 
rgr_cmd_15spp_allses <- bind_cols(rgr_cmd_15spp_means, rgr_cmd_15spp_lowCI, rgr_cmd_15spp_upperCI)
#select subset of columns and rename because the naming changed during merging 
rgr_cmd_15spp_allses <- dplyr::select(rgr_cmd_15spp_allses, spp...1, spp.1...2, predictor...3, id...5, mean, low_CI, upper_CI) %>% 
  rename(rep = spp...1, spp = spp.1...2, id = id...5, predictor = predictor...3)


#reorder the species according to relative frequency 
#this doesn't objectively look like it changes anything, but it does plot in the right order
rgr_cmd_15spp_allses$spp <- factor(rgr_cmd_15spp_allses$spp, levels=rev(c("Pseudotsuga menziesii", "Pinus ponderosa",
                                                                            "Tsuga heterophylla", "Pinus contorta", "Abies concolor", "Abies amabilis",
                                                                            "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana",
                                                                            "Thuja plicata", "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis",
                                                                            "Calocedrus decurrens", "Juniperus occidentalis"))) # plot species by rel freq order

### plot the bootstrapped coefficients per species
rgr_cmd_15 <- ggplot(data = rgr_cmd_15spp_allses, aes(y = spp, x = mean, xmin=low_CI, xmax=upper_CI)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.5, color= "black", lwd=0.5) +
  geom_vline(xintercept=0, size=0.5, linetype='dashed') +
  facet_grid( ~ predictor, scales = 'fixed') +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-9e-04, 6e-04) + #put on same scale as the fire plot 
  ggtitle("") +
  theme_bw()

rgr_cmd_15 + theme(axis.text.y = element_text(face="italic", size = 12)) + theme(axis.text.x = element_text(size = 12)) + theme(strip.text.x = element_text(size = 12))


