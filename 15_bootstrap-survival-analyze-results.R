#############################################################################################################
#
#   Tree Responses to Regional Gradients -- Bootstrapping survival mixed-effects models & analysis of results 
#
#   L. McKinley Nevins, laura.nevins@wsu.edu, 26 June 2023
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

### load data
setwd("~/Dropbox/WSU/Cascade_Project/Kamiak/Resub2")
survival <- read.csv("survival_full.csv")

#################################################################################################################
##Generalized Linear Mixed Models 

#perform once to get coefficients 

#fire#####
start_time <- Sys.time()
m_surv_fire <- lme4::glmer(surv ~ firelog * pc1            # second-order interaction term
                           + (1 + firelog * pc1 | spp)      # separate slopes for species 
                           + (dia_y1 + I(dia_y1^2))               # standardized tree size 
                           + (1 | census_int)                           # interval between censuses 
                           + (1 | censusyr)                       # separate means for year
                           + (1 | plt_cn),                        # separate means for plot
                           data      = survival,  
                           binomial(link="logit"),
                           na.action = 'na.fail',                 # NA handling
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e9)))                 
end_time <- Sys.time()
end_time - start_time # time elapsed ~ 20.25 minutes for standard binomial survival 

#check for singularity
isSingular(m_surv_fire, tol = 1e-4)


#summaries
summary(m_surv_fire)
coef(m_surv_fire)
car::Anova(m_surv_fire, type=3)
performance::r2_nakagawa(m_surv_fire,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)


##cmd#####
start_time <- Sys.time()
m_surv_cmd <- lme4::glmer(surv ~ cmdlog * pc1             # second-order interaction terms
                          + (1 + cmdlog * pc1 | spp)      # separate slopes for species 
                          + (dia_y1 + I(dia_y1^2))               # standardized tree size 
                          + (1 | census_int)                         # interval between censuses 
                          + (1 | censusyr)                       # separate means for year
                          + (1 | plt_cn),                        # separate means for plot
                          data      = survival, 
                          binomial(link="logit"),
                          na.action = 'na.fail',                   # NA handling
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e8)))    
end_time <- Sys.time()
end_time - start_time


#check for singularity
isSingular(m_surv_cmd, tol = 1e-4)


#summaries
summary(m_surv_cmd)
coef(m_surv_cmd)
car::Anova(m_surv_cmd, type=3)
performance::r2_nakagawa(m_surv_cmd,
                         by_group   = FALSE,
                         tolerance  = 1e-05,
                         ci         = NULL,
                         iterations = 99)


######### Bootstrapping ############################################################

##Survival - Fire x PC1##################

#takes ~24 days to run 
nboots <- 999
cxs <- coef(m_surv_fire)$spp[,c('firelog','pc1','firelog:pc1', 'I(dia_y1^2)')] # per-species coefficients for two main effects and interaction
a <- array(NA, dim=c(dim(cxs),nboots), dimnames=list(rownames(cxs),colnames(cxs),NULL))
# ! ! ! ! TIMEWARN ! ! ! ! about 20 min per iteration ! ! ! ! ! ! ! ! ! ! ! ! ! !
for (b in 1:nboots) {
  cat('iteration', b, 'of', nboots, '....... \n')
  nr <- NROW(survival)
  dx <- survival[sample(1:nr, size=nr, replace=FALSE),,] # draw with replacement
  m_surv_fire <- lme4::glmer(surv ~ firelog * pc1            # second-order interaction term
                             + (1 + firelog * pc1 | spp)      # separate slopes for species 
                             + (dia_y1 + I(dia_y1^2))               # standardized tree size 
                             + (1 | census_int)                         # interval between censuses 
                             + (1 | censusyr)                       # separate means for year
                             + (1 | plt_cn),                        # separate means for plot
                             data      = dx,  
                             binomial(link="logit"),
                             na.action = 'na.fail',                 # NA handling
                             control=glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun=2e5)))    
  a[,,b] <- as.matrix(coef(m_surv_fire)$spp[,c('firelog','pc1','firelog:pc1', 'I(dia_y1^2)')])
}


#move the results array into a clearer name
surv_fire_mod999 <- a

#save the results
setwd("~/Dropbox/WSU/Cascade_Project/J of Ecol Submission/Resubmission2")
save(surv_fire_mod999, file="surv_fire_mod999.rda")



##Survival - CMD x PC1##################

#takes ~30 days to run 
nboots <- 999
cxs <- coef(m_surv_cmd)$spp[,c('cmdlog','pc1','cmdlog:pc1','I(dia_y1^2)')] # per-species coefficients for two main effects, interaction, and diameter
a <- array(NA, dim=c(dim(cxs),nboots), dimnames=list(rownames(cxs),colnames(cxs),NULL))
# ! ! ! ! TIMEWARN ! ! ! ! about 30 min per iteration ! ! ! ! ! ! ! ! ! ! ! ! ! !
for (b in 1:nboots) {
  cat('iteration', b, 'of', nboots, '....... \n')
  nr <- NROW(survival)
  dx <- survival[sample(1:nr, size=nr, replace=FALSE),,] # draw with replacement
  start_time <- Sys.time()
  m_surv_cmd <- lme4::glmer(surv ~ cmdlog * pc1             # second-order interaction terms
                            + (1 + cmdlog * pc1 | spp)      # separate slopes for species 
                            + (dia_y1 + I(dia_y1^2))               # standardized tree size 
                            + (1 | census_int)                           # interval between censuses 
                            + (1 | censusyr)                       # separate means for year
                            + (1 | plt_cn),                        # separate means for plot
                            data      = dx, 
                            binomial(link="logit"),
                            na.action = 'na.fail',                   # NA handling
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e7)))    
  end_time <- Sys.time()
  end_time - start_time
  a[,,b] <- as.matrix(coef(m_surv_cmd)$spp[,c('cmdlog','pc1','cmdlog:pc1','I(dia_y1^2)')])
}

#move the results array into a clearer name
surv_cmd_mod999 <- a

#save the results
setwd("~/Dropbox/WSU/Cascade_Project/J of Ecol Submission/Resubmission2")
save(surv_cmd_mod999, file="surv_cmd_mod999.rda")

####################################################################################################
##Load bootstrapped results arrays 

setwd("~/Documents/McKinley/Cascade_Project/")

#fire model 
load(file = 'surv_fire_mod999.rda')

#cmd model 
load(file = 'surv_cmd_mod999.rda')

#########Generating Final Plots from the Results###################################################

####### Fire Results###########

### reshape to get mean coefficients
surv_fire_allspp <- surv_fire_mod999          # resulting coefficients array in full
# reshape it
surv_fire_allspp <- matrix(aperm(surv_fire_allspp, c(1, 3, 2)), nrow = dim(surv_fire_allspp)[1] * dim(surv_fire_allspp)[3], dimnames = list(rep(rownames(surv_fire_allspp),dim(surv_fire_allspp)[3]), colnames(surv_fire_allspp)))
surv_fire_allspp <- data.frame(spp = rownames(surv_fire_allspp), surv_fire_allspp, row.names=NULL)
surv_fire_allspp <- surv_fire_allspp[order(surv_fire_allspp$spp),,]

#55,944 lines in the datafile, because there were 56 species iterated 999 times 

#get mean, sd, and se for values across all species
surv_fire_allspp_vals <- surv_fire_allspp %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, se = ~sd(.x)/sqrt(length(.x)))))

####Calculate confidence intervals  

#for fire
firehist <- hist(surv_fire_allspp$firelog)
quantile(surv_fire_allspp$firelog, c(0.025, 0.975))

#set lower and upper bound from output 
firelog_low_CI <- -0.2597152
firelog_upper_CI <- 0.2817277

#as dataframes
firelog_low_CI <- as.data.frame(firelog_low_CI)
firelog_upper_CI <- as.data.frame(firelog_upper_CI)


##for pc1
pc1hist <- hist(surv_fire_allspp$pc1)
quantile(surv_fire_allspp$pc1, c(0.025, 0.975))

# set lower bound and upper bound
pc1_low_CI <- -0.1160388
pc1_upper_CI <- 0.3452171

#as dataframes
pc1_low_CI <- as.data.frame(pc1_low_CI)
pc1_upper_CI <- as.data.frame(pc1_upper_CI)


##for firelog x pc1
fire.pc1hist <- hist(surv_fire_allspp$firelog.pc1)
quantile(surv_fire_allspp$firelog.pc1, c(0.025, 0.975))

#set lower bound and upper bound
firelog.pc1_low_CI <- -0.1933320
firelog.pc1_upper_CI <- 0.1869505

#as dataframes
firelog.pc1_low_CI <- as.data.frame(firelog.pc1_low_CI)
firelog.pc1_upper_CI <- as.data.frame(firelog.pc1_upper_CI)

##for diameter
diameterhist <- hist(surv_fire_allspp$I.dia_y1.2.)
quantile(surv_fire_allspp$I.dia_y1.2., c(0.025, 0.975))

# Calculating lower bound and upper bound
dia_low_CI <- -0.1497541
dia_upper_CI <- -0.1497541

#as dataframes
dia_low_CI <- as.data.frame(dia_low_CI)
dia_upper_CI <- as.data.frame(dia_upper_CI)


#subset surv_fire_allspp_vals to get the means 
surv_fire_allspp_vals <- surv_fire_allspp_vals %>%
  dplyr::select(firelog_mean, pc1_mean, firelog.pc1_mean, I.dia_y1.2._mean)


#bind everything together 
allspp_surv_fire_intervals <- bind_cols(surv_fire_allspp_vals, firelog_low_CI, firelog_upper_CI, 
                                       pc1_low_CI, pc1_upper_CI,
                                       firelog.pc1_low_CI, firelog.pc1_upper_CI,
                                       dia_low_CI, dia_upper_CI)

#combine into a new dataframe 
surv_fire_allspp_vals <- data.frame(coefficient = c('Fire', 'Funct Neigh', 'Fire x Funct Neigh', 'Initial Tree Size'),
                                   mean = c(0.0317073, 0.09226432, 0.01431195, -0.1496915),
                                   low_confint = c(-0.2597152, -0.1160388, -0.193332, -0.1496915),
                                   upper_confint = c(0.2817277, 0.3452171, 0.1869505, -0.1496915),
                                   sig = c('no', 'no', 'no', 'yes'))


#plot coefficients and confidence intervals for all species combined 
allplot_survfire <-ggplot(surv_fire_allspp_vals,aes(x=mean, y=coefficient)) +
  geom_point() + 
  geom_errorbar(aes(xmin=low_confint, xmax=upper_confint), colour="black", width=.1) +
  geom_point(color='black', shape=21, size=2, aes(fill=factor(sig)), show.legend = FALSE) + 
  scale_fill_manual(values=c('white', 'black')) +
  geom_line() +
  geom_vline(xintercept=0, size=1) +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.4, 0.5) +
  ggtitle("") +
  theme_bw()


allplot_survfire + theme(axis.text.y = element_text(size = 10), 
                        axis.text.x = element_text(size = 10),)


######## assess each species individually #####################################################################################

### reshape the resulting bootstrapped coefficients array
surv_fire_allspp <- surv_fire_mod999                  # resulting coefficients array in full
# reshape it
surv_fire_15spp <- matrix(aperm(surv_fire_allspp, c(1, 3, 2)), nrow = dim(surv_fire_allspp)[1] * dim(surv_fire_allspp)[3], dimnames = list(rep(rownames(surv_fire_allspp),dim(surv_fire_allspp)[3]), colnames(surv_fire_allspp)))
surv_fire_15spp <- data.frame(spp = rownames(surv_fire_15spp), surv_fire_15spp, row.names=NULL)
surv_fire_15spp <- surv_fire_15spp[order(surv_fire_15spp$spp),,]

#subset 15 most abundant species 
surv_fire_15spp = surv_fire_15spp[surv_fire_15spp$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                                          "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                                          "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 
#reduces to 14,985 observations for the 5 variables 

#get mean for each species 
surv_fire_stats <- surv_fire_15spp %>% 
  group_by(spp) %>% 
  summarise(across(everything(), list(mean = mean)))

#check the distributions of the coefficient values 
hist(surv_fire_15spp$firelog)
hist(surv_fire_15spp$pc1)
hist(surv_fire_15spp$firelog.pc1)

#this calculates the upper and lower confidence intervals (here as the 2.5% and 97.5% quantiles, and organizes them by species)

#for fire
firelog_quants <- surv_fire_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(firelog, c(0.025, 0.975)), "quantile", "firelog"))

#for pc1
pc1_quants <- surv_fire_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(pc1, c(0.025, 0.975)), "quantile", "pc1"))

#for fire x pc1
firelog.pc1_quants <- surv_fire_15spp %>%
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
surv_fire_15spp_means <- dplyr::select(.data = surv_fire_stats, spp, firelog_mean, pc1_mean, firelog.pc1_mean)

surv_fire_15spp_means <- surv_fire_15spp_means %>% 
  rename(Fire = firelog_mean, Funct_Neigh = pc1_mean, FirexFunct_Neigh = firelog.pc1_mean)

#reshape the means dataset 
surv_fire_15spp_means <- data.frame(spp = rownames(surv_fire_15spp_means), surv_fire_15spp_means, row.names=NULL)
surv_fire_15spp_means <- surv_fire_15spp_means[order(surv_fire_15spp_means$spp),,]
surv_fire_15spp_means  <- reshape(surv_fire_15spp_means, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='mean',
                                 times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(surv_fire_15spp_means)))

#reshape the low CI dataset 
surv_fire_15spp_lowCI <- data.frame(spp = rownames(low_quants), low_quants, row.names=NULL)
surv_fire_15spp_lowCI <- surv_fire_15spp_lowCI[order(surv_fire_15spp_lowCI$spp),,]
surv_fire_15spp_lowCI  <- reshape(surv_fire_15spp_lowCI, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='low_CI',
                                 times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                 new.row.names=1:prod(dim(surv_fire_15spp_lowCI)))

#reshape the upper CI dataset 
surv_fire_15spp_upperCI <- data.frame(spp = rownames(upper_quants), upper_quants, row.names=NULL)
surv_fire_15spp_upperCI <- surv_fire_15spp_upperCI[order(surv_fire_15spp_upperCI$spp),,]
surv_fire_15spp_upperCI  <- reshape(surv_fire_15spp_upperCI, dir='long', varying=c('Fire','Funct_Neigh','FirexFunct_Neigh'), v.names ='upper_CI',
                                   times=c('Fire','Funct_Neigh','FirexFunct_Neigh'), timevar='predictor',
                                   new.row.names=1:prod(dim(surv_fire_15spp_upperCI)))

#merge dataframes to get means and CI's all together 
surv_fire_15spp_allses <- bind_cols(surv_fire_15spp_means, surv_fire_15spp_lowCI, surv_fire_15spp_upperCI)
#select subset of columns and rename because the naming changed during merging 
#there's probably a cleaner way to do this but I don't know it
surv_fire_15spp_allses <- dplyr::select(surv_fire_15spp_allses, spp...1, spp.1...2, predictor...3, id...5, mean, low_CI, upper_CI) %>% 
  rename(rep = spp...1, spp = spp.1...2, id = id...5, predictor = predictor...3)

#reorder the species according to relative frequency 
#this doesn't objectively look like it changes anything, but it does plot in the right order
surv_fire_15spp_allses$spp <- factor(surv_fire_15spp_allses$spp, levels=rev(c("Pseudotsuga menziesii", "Pinus ponderosa",
                                                                            "Tsuga heterophylla", "Pinus contorta", "Abies concolor", "Abies amabilis",
                                                                            "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana",
                                                                            "Thuja plicata", "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis",
                                                                            "Calocedrus decurrens", "Juniperus occidentalis"))) # plot species by rel freq order

#don't need to change the color of any of these points because all of the species are significant, 
#and none cross the zero line 

### plot the bootstrapped coefficients per species

surv_fire_15 <- ggplot(data = surv_fire_15spp_allses, aes(y = spp, x = mean, xmin=low_CI, xmax=upper_CI)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.5, color= "black", lwd=0.5) +
  geom_vline(xintercept=0, size=0.5, linetype='dashed') +
  facet_grid( ~ predictor, scales = 'fixed') +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.4, 0.6) +
  ggtitle("") +
  theme_bw()

surv_fire_15 + theme(axis.text.y = element_text(face="italic", size = 10)) + theme(axis.text.x = element_text(size = 10))


###########################################################################################################################

####### CMD Results###########

### reshape to get mean coefficients for each species 
surv_cmd_allspp <- surv_cmd_mod999          # resulting coefficients array in full
# reshape it
surv_cmd_allspp <- matrix(aperm(surv_cmd_allspp, c(1, 3, 2)), nrow = dim(surv_cmd_allspp)[1] * dim(surv_cmd_allspp)[3], dimnames = list(rep(rownames(surv_cmd_allspp),dim(surv_cmd_allspp)[3]), colnames(surv_cmd_allspp)))
surv_cmd_allspp <- data.frame(spp = rownames(surv_cmd_allspp), surv_cmd_allspp, row.names=NULL)
surv_cmd_allspp <- surv_cmd_allspp[order(surv_cmd_allspp$spp),,]

#55,944 lines in the datafile, because there were 56 species iterated 999 times 

#get mean, sd, and se for values across all species
surv_cmd_allspp_vals <- surv_cmd_allspp %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, se = ~sd(.x)/sqrt(length(.x)))))

####Calculating confidence intervals

#for cmd
cmdhist <- hist(surv_cmd_allspp$cmdlog)
quantile(surv_cmd_allspp$cmdlog, c(0.025, 0.975))

#set lower and upper bound from output 
cmdlog_low_CI <- -0.3331596
cmdlog_upper_CI <- 0.4814786

#as dataframes
cmdlog_low_CI <- as.data.frame(cmdlog_low_CI)
cmdlog_upper_CI <- as.data.frame(cmdlog_upper_CI)


##for pc1
pc1.2hist <- hist(surv_cmd_allspp$pc1)
quantile(surv_cmd_allspp$pc1, c(0.025, 0.975))

# set lower bound and upper bound
pc1_low_CI <- -0.1236591
pc1_upper_CI <- 0.2863100

#as dataframes
pc1_low_CI <- as.data.frame(pc1_low_CI)
pc1_upper_CI <- as.data.frame(pc1_upper_CI)


##for cmdlog x pc1
cmd.pc1hist <- hist(surv_cmd_allspp$cmdlog.pc1)
quantile(surv_cmd_allspp$cmdlog.pc1, c(0.025, 0.975))

#set lower bound and upper bound
cmdlog.pc1_low_CI <- -0.04934949
cmdlog.pc1_upper_CI <- 0.22824866

#as dataframes
cmdlog.pc1_low_CI <- as.data.frame(cmdlog.pc1_low_CI)
cmdlog.pc1_upper_CI <- as.data.frame(cmdlog.pc1_upper_CI)


##for diameter
diameterhist <- hist(surv_cmd_allspp$I.dia_y1.2.)
quantile(surv_cmd_allspp$I.dia_y1.2., c(0.025, 0.975))

# Calculating lower bound and upper bound
dia_low_CI <- -0.1509592
dia_upper_CI <- -0.1508443

#as dataframes
dia_low_CI <- as.data.frame(dia_low_CI)
dia_upper_CI <- as.data.frame(dia_upper_CI)


#subset surv_cmd_allspp_vals to get the means 
surv_cmd_allspp_vals <- surv_cmd_allspp_vals %>%
  dplyr::select(cmdlog_mean, pc1_mean, cmdlog.pc1_mean, I.dia_y1.2._mean)


#bind everything together 
allspp_surv_cmd_intervals <- bind_cols(surv_cmd_allspp_vals, cmdlog_low_CI, cmdlog_upper_CI, 
                                      pc1_low_CI, pc1_upper_CI,
                                      cmdlog.pc1_low_CI, cmdlog.pc1_upper_CI,
                                      dia_low_CI, dia_upper_CI)

#combine into a new dataframe 
surv_cmd_allspp_vals <- data.frame(coefficient = c('CMD', 'Funct Neigh', 'CMD x Funct Neigh', 'Initial Tree Size'),
                                  mean = c(0.09120303, 0.05780241, 0.09617802, -0.1508584),
                                  low_confint = c(-0.3331596, -0.1236591, -0.04934949, -0.1509592),
                                  upper_confint = c(0.4814786, 0.28631, 0.2282487, -0.1508443),
                                  sig = c('no', 'no', 'no', 'yes'))


#plot coefficients and confidence intervals for all species combined 
allplot_survcmd <-ggplot(surv_cmd_allspp_vals,aes(x=mean,y=coefficient)) +
  geom_point() + 
  geom_errorbar(aes(xmin=low_confint, xmax=upper_confint), colour="black", width=.1) +
  geom_point(color='black', shape=21, size=2, aes(fill=factor(sig)), show.legend = FALSE) + 
  scale_fill_manual(values=c('white', 'black')) +
  geom_line() +
  geom_vline(xintercept=0, size=1) +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.4, 0.5) + #put on same scale as the fire plot 
  ggtitle("") +
  theme_bw()


allplot_survcmd + theme(axis.text.y = element_text(size = 10), 
                       axis.text.x = element_text(size = 10),)

######## assess each species individually #####################################################################################

### reshape the resulting bootstrapped coefficients array
surv_cmd_allspp <- surv_cmd_mod999                  # resulting coefficients array in full
# reshape it
surv_cmd_15spp <- matrix(aperm(surv_cmd_allspp, c(1, 3, 2)), nrow = dim(surv_cmd_allspp)[1] * dim(surv_cmd_allspp)[3], dimnames = list(rep(rownames(surv_cmd_allspp),dim(surv_cmd_allspp)[3]), colnames(surv_cmd_allspp)))
surv_cmd_15spp <- data.frame(spp = rownames(surv_cmd_15spp), surv_cmd_15spp, row.names=NULL)
surv_cmd_15spp <- surv_cmd_15spp[order(surv_cmd_15spp$spp),,]

#subset 15 most abundant species 
surv_cmd_15spp = surv_cmd_15spp[surv_cmd_15spp$spp %in% c("Pseudotsuga menziesii", "Pinus ponderosa", "Tsuga heterophylla", "Pinus contorta", "Abies concolor", 
                                                       "Abies amabilis", "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana", "Thuja plicata", 
                                                       "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis", "Calocedrus decurrens", "Juniperus occidentalis"), ] 
#reduces to 14,985 observations for the 5 variables 

#get mean for each species 
surv_cmd_stats <- surv_cmd_15spp %>% 
  group_by(spp) %>% 
  summarise(across(everything(), list(mean = mean)))

#check the distributions of the coefficient values 
hist(surv_cmd_15spp$cmdlog)
hist(surv_cmd_15spp$pc1)
hist(surv_cmd_15spp$cmdlog.pc1)

#this calculates the upper and lower confidence intervals (here as the 2.5% and 97.5% quantiles, and organizes them by species)

#for cmd
cmdlog_quants <- surv_cmd_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(cmdlog, c(0.025, 0.975)), "quantile", "cmdlog"))

#for pc1
pc1_quants <- surv_cmd_15spp %>%
  group_by(spp) %>%
  summarise(enframe(quantile(pc1, c(0.025, 0.975)), "quantile", "pc1"))

#for cmd x pc1
cmdlog.pc1_quants <- surv_cmd_15spp %>%
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
surv_cmd_15spp_means <- dplyr::select(.data = surv_cmd_stats, spp, cmdlog_mean, pc1_mean, cmdlog.pc1_mean)

surv_cmd_15spp_means <- surv_cmd_15spp_means %>% 
  rename(CMD = cmdlog_mean, Funct_Neigh = pc1_mean, CMDxFunct_Neigh = cmdlog.pc1_mean)

#reshape the means dataset 
surv_cmd_15spp_means <- data.frame(spp = rownames(surv_cmd_15spp_means), surv_cmd_15spp_means, row.names=NULL)
surv_cmd_15spp_means <- surv_cmd_15spp_means[order(surv_cmd_15spp_means$spp),,]
surv_cmd_15spp_means  <- reshape(surv_cmd_15spp_means, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='mean',
                                times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                new.row.names=1:prod(dim(surv_cmd_15spp_means)))

#reshape the low CI dataset 
surv_cmd_15spp_lowCI <- data.frame(spp = rownames(low_quants), low_quants, row.names=NULL)
surv_cmd_15spp_lowCI <- surv_cmd_15spp_lowCI[order(surv_cmd_15spp_lowCI$spp),,]
surv_cmd_15spp_lowCI  <- reshape(surv_cmd_15spp_lowCI, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='low_CI',
                                times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                new.row.names=1:prod(dim(surv_cmd_15spp_lowCI)))

#reshape the upper CI dataset 
surv_cmd_15spp_upperCI <- data.frame(spp = rownames(upper_quants), upper_quants, row.names=NULL)
surv_cmd_15spp_upperCI <- surv_cmd_15spp_upperCI[order(surv_cmd_15spp_upperCI$spp),,]
surv_cmd_15spp_upperCI  <- reshape(surv_cmd_15spp_upperCI, dir='long', varying=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), v.names ='upper_CI',
                                  times=c('CMD','Funct_Neigh','CMDxFunct_Neigh'), timevar='predictor',
                                  new.row.names=1:prod(dim(surv_cmd_15spp_upperCI)))


#merge dataframes to get means and CI's all together 
surv_cmd_15spp_allses <- bind_cols(surv_cmd_15spp_means, surv_cmd_15spp_lowCI, surv_cmd_15spp_upperCI)
#select subset of columns and rename because the naming changed during merging 
surv_cmd_15spp_allses <- dplyr::select(surv_cmd_15spp_allses, spp...1, spp.1...2, predictor...3, id...5, mean, low_CI, upper_CI) %>% 
  rename(rep = spp...1, spp = spp.1...2, id = id...5, predictor = predictor...3)


#reorder the species according to relative frequency 
#this doesn't objectively look like it changes anything, but it does plot in the right order
surv_cmd_15spp_allses$spp <- factor(surv_cmd_15spp_allses$spp, levels=rev(c("Pseudotsuga menziesii", "Pinus ponderosa",
                                                                          "Tsuga heterophylla", "Pinus contorta", "Abies concolor", "Abies amabilis",
                                                                          "Abies grandis", "Lithocarpus densiflorus", "Tsuga mertensiana",
                                                                          "Thuja plicata", "Alnus rubra", "Abies lasiocarpa", "Quercus chrysolepis",
                                                                          "Calocedrus decurrens", "Juniperus occidentalis"))) # plot species by rel freq order

### plot the bootstrapped coefficients per species
surv_cmd_15 <- ggplot(data = surv_cmd_15spp_allses, aes(y = spp, x = mean, xmin=low_CI, xmax=upper_CI)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.5, color= "black", lwd=0.5) +
  geom_vline(xintercept=0, size=0.5, linetype='dashed') +
  facet_grid( ~ predictor, scales = 'fixed') +
  labs(x = "Effect Size (Mean Bootstrapped Coefficient Value)", y = "") + xlim(-0.4, 0.6) + #put on same scale as the fire plot 
  ggtitle("") +
  theme_bw()

surv_cmd_15 + theme(axis.text.y = element_text(face="italic", size = 10)) + theme(axis.text.x = element_text(size = 10))

