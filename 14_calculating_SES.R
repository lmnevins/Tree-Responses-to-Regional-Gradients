library(data.table)
library(purrr)
library(dplyr)

#Load barck thickness randomization data
setwd("~/Dropbox/Functional_neighborhood_randomizations")

temp <- list.files(pattern='functi_neigh')

bark <- map(temp, fread, select = 'ncis_Bark_thickness') %>%
  bind_cols

#Load observed bark thickness data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
bark_obs <- neigh_obs[,c("focal","ncis_Bark_thickness")]

#Add focal to the bark data
focal <- neigh_obs[,1]
bark_focal <- bind_cols(focal, bark)

#Calculate mean and SD from the bark data
mean_bark <- sapply(1:nrow(bark_focal), FUN = function(i){mean(as.numeric(bark_focal[i,-1]), na.rm = T)})
sd_bark <- sapply(1:nrow(bark_focal), FUN = function(i){sd(as.numeric(bark_focal[i,-1]), na.rm = T)})
bark_mean_sd <- cbind.data.frame(bark_focal[,1], mean_bark, sd_bark)
colnames(bark_mean_sd)[1] <- "focal"


bark_ses <- inner_join(bark_obs, bark_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_Bark_thickness - mean_bark) / sd_bark ) %>%
  select(focal, stand_size_bark = stand_size)

setwd("~/Dropbox/Cascade_project/Results/Functional_neighborhood_SES")

write.csv(bark_ses, "ses_bark_thickness.csv", row.names = FALSE)

##############################LITTER DECOMPOSITION##########################
litter <- map(temp, fread, select = 'ncis_Litter_decomposition') %>%
  bind_cols

#Load observed litter data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
litter_obs <- neigh_obs[,c("focal","ncis_Litter_decomposition")]

#Add focal to the litter data
focal <- neigh_obs[,1]
litter_focal <- bind_cols(focal, litter)

#Calculate mean and SD from the litter data
mean_litter <- sapply(1:nrow(litter_focal), FUN = function(i){mean(as.numeric(litter_focal[i,-1]), na.rm = T)})
sd_litter <- sapply(1:nrow(litter_focal), FUN = function(i){sd(as.numeric(litter_focal[i,-1]), na.rm = T)})
litter_mean_sd <- cbind.data.frame(litter_focal[,1], mean_litter, sd_litter)
colnames(litter_mean_sd)[1] <- "focal"

litter_ses <- inner_join(litter_obs, litter_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_Litter_decomposition - mean_litter) / sd_litter ) %>%
  select(focal, stand_size_litter = stand_size)

write.csv(litter_ses, "ses_litter_decomposition.csv", row.names = FALSE)


##############################P50N##########################
P50 <- map(temp, fread, select = 'ncis_P50') %>%
  bind_cols

#Load observed P50 data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
P50_obs <- neigh_obs[,c("focal","ncis_P50")]

#Add focal to the P50 data
focal <- neigh_obs[,1]
P50_focal <- bind_cols(focal, P50)

#Calculate mean and SD from the P50 data
mean_P50 <- sapply(1:nrow(P50_focal), FUN = function(i){mean(as.numeric(P50_focal[i,-1]), na.rm = T)})
sd_P50 <- sapply(1:nrow(P50_focal), FUN = function(i){sd(as.numeric(P50_focal[i,-1]), na.rm = T)})
P50_mean_sd <- cbind.data.frame(P50_focal[,1], mean_P50, sd_P50)
colnames(P50_mean_sd)[1] <- "focal"

P50_ses <- inner_join(P50_obs, P50_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_P50 - mean_P50) / sd_P50) %>%
  select(focal, stand_size_P50 = stand_size)

write.csv(P50_ses, "ses_P50.csv", row.names = FALSE)

##############################STEM CONDUCTIVITY##########################
conductivity <- map(temp, fread, select = 'ncis_Stem_conductivity') %>%
  bind_cols

#Load observed conductivity data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
conductivity_obs <- neigh_obs[,c("focal","ncis_Stem_conductivity")]

#Add focal to the conductivity data
focal <- neigh_obs[,1]
conductivity_focal <- bind_cols(focal, conductivity)

#Calculate mean and SD from the conductivity data
mean_conductivity <- sapply(1:nrow(conductivity_focal), FUN = function(i){mean(as.numeric(conductivity_focal[i,-1]), na.rm = T)})
sd_conductivity <- sapply(1:nrow(conductivity_focal), FUN = function(i){sd(as.numeric(conductivity_focal[i,-1]), na.rm = T)})
conductivity_mean_sd <- cbind.data.frame(conductivity_focal[,1], mean_P50, sd_P50)
colnames(conductivity_mean_sd)[1] <- "focal"


conductivity_ses <- inner_join(conductivity_obs, conductivity_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_Stem_conductivity - mean_conductivity) / sd_conductivity) %>%
  select(focal, stand_size_conductivity = stand_size)

write.csv(conductivity_ses, "ses_Stem_conductivity.csv", row.names = FALSE)

##############################ROOT##########################
root <- map(temp, fread, select = 'ncis_Root_depth') %>%
  bind_cols

#Load observed root data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
root_obs <- neigh_obs[,c("focal","ncis_Root_depth")]

#Add focal to the root data
focal <- neigh_obs[,1]
root_focal <- bind_cols(focal, root)

#Calculate mean and SD from the conductivity data
mean_root <- sapply(1:nrow(root_focal), FUN = function(i){mean(as.numeric(root_focal[i,-1]), na.rm = T)})
sd_root <- sapply(1:nrow(root_focal), FUN = function(i){sd(as.numeric(root_focal[i,-1]), na.rm = T)})
root_mean_sd <- cbind.data.frame(root_focal[,1], mean_root, sd_root)
colnames(root_mean_sd)[1] <- "focal"

root_ses <- inner_join(root_obs, root_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_Root_depth - mean_root) / sd_root) %>%
  select(focal, stand_size_root = stand_size)

write.csv(root_ses, "ses_Root_depth.csv", row.names = FALSE)

##############################PC1 FIRE##########################
pc1_fire <- map(temp, fread, select = 'ncis_PC1_fire') %>%
  bind_cols

#Load observed fire data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
pc1_fire_obs <- neigh_obs[,c("focal","ncis_PC1_fire")]

#Add focal to the fire data
focal <- neigh_obs[,1]
pc1_fire_focal <- bind_cols(focal, pc1_fire)

#Calculate mean and SD from the conductivity data
mean_pc1_fire <- sapply(1:nrow(pc1_fire_focal), FUN = function(i){mean(as.numeric(pc1_fire_focal[i,-1]), na.rm = T)})
sd_pc1_fire <- sapply(1:nrow(pc1_fire_focal), FUN = function(i){sd(as.numeric(pc1_fire_focal[i,-1]), na.rm = T)})
pc1_fire_mean_sd <- cbind.data.frame(pc1_fire_focal[,1], mean_pc1_fire, sd_pc1_fire)
colnames(pc1_fire_mean_sd)[1] <- "focal"

pc1_fire_ses <- inner_join(pc1_fire_obs, pc1_fire_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_PC1_fire - mean_pc1_fire) / sd_pc1_fire) %>%
  select(focal, stand_size_pc1_fire = stand_size)

write.csv(pc1_fire_ses, "ses_PC1_fire.csv", row.names = FALSE)

##############################PC2 FIRE##########################
pc2_fire <- map(temp, fread, select = 'ncis_PC2_fire') %>%
  bind_cols

#Load observed fire data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
pc2_fire_obs <- neigh_obs[,c("focal","ncis_PC2_fire")]

#Add focal to the fire data
focal <- neigh_obs[,1]
pc2_fire_focal <- bind_cols(focal, pc2_fire)

#Calculate mean and SD from the fire data
mean_pc2_fire <- sapply(1:nrow(pc2_fire_focal), FUN = function(i){mean(as.numeric(pc2_fire_focal[i,-1]), na.rm = T)})
sd_pc2_fire <- sapply(1:nrow(pc2_fire_focal), FUN = function(i){sd(as.numeric(pc2_fire_focal[i,-1]), na.rm = T)})
pc2_fire_mean_sd <- cbind.data.frame(pc2_fire_focal[,1], mean_pc2_fire, sd_pc2_fire)
colnames(pc2_fire_mean_sd)[1] <- "focal"

pc2_fire_ses <- inner_join(pc2_fire_obs, pc2_fire_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_PC2_fire - mean_pc2_fire) / sd_pc2_fire) %>%
  select(focal, stand_size_pc2_fire = stand_size)

write.csv(pc2_fire_ses, "ses_PC2_fire.csv", row.names = FALSE)

##############################PC1 DROUGHT##########################
pc1_drought <- map(temp, fread, select = 'ncis_PC1_drought') %>%
  bind_cols

#Load observed drought data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
pc1_drought_obs <- neigh_obs[,c("focal","ncis_PC1_drought")]

#Add focal to the drought data
focal <- neigh_obs[,1]
pc1_drought_focal <- bind_cols(focal, pc1_drought)

#Calculate mean and SD from the drought data
mean_pc1_drought <- sapply(1:nrow(pc1_drought_focal), FUN = function(i){mean(as.numeric(pc1_drought_focal[i,-1]), na.rm = T)})
sd_pc1_drought <- sapply(1:nrow(pc1_drought_focal), FUN = function(i){sd(as.numeric(pc1_drought_focal[i,-1]), na.rm = T)})
pc1_drought_mean_sd <- cbind.data.frame(pc1_drought_focal[,1], mean_pc1_drought, sd_pc1_drought)
colnames(pc1_drought_mean_sd)[1] <- "focal"

pc1_drought_ses <- inner_join(pc1_drought_obs, pc1_drought_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_PC1_drought - mean_pc1_drought) / sd_pc1_drought) %>%
  select(focal, stand_size_pc1_drought = stand_size)

write.csv(pc1_drought_ses, "ses_PC1_drought.csv", row.names = FALSE)

##############################PC2 DROUGHT##########################
pc2_drought <- map(temp, fread, select = 'ncis_PC2_drought') %>%
  bind_cols

#Load observed drought data
neigh_obs <- fread("functi_neigh_obs.csv", header=TRUE)
pc2_drought_obs <- neigh_obs[,c("focal","ncis_PC2_drought")]

#Add focal to the drought data
focal <- neigh_obs[,1]
pc2_drought_focal <- bind_cols(focal, pc2_drought)

#Calculate mean and SD from the drought data
mean_pc2_drought <- sapply(1:nrow(pc2_drought_focal), FUN = function(i){mean(as.numeric(pc2_drought_focal[i,-1]), na.rm = T)})
sd_pc2_drought <- sapply(1:nrow(pc2_drought_focal), FUN = function(i){sd(as.numeric(pc2_drought_focal[i,-1]), na.rm = T)})
pc2_drought_mean_sd <- cbind.data.frame(pc2_drought_focal[,1], mean_pc2_drought, sd_pc2_drought)
colnames(pc2_drought_mean_sd)[1] <- "focal"

pc2_drought_ses <- inner_join(pc2_drought_obs, pc2_drought_mean_sd, by = "focal") %>%
  mutate(stand_size = (ncis_PC2_drought - mean_pc2_drought) / sd_pc2_drought) %>%
  select(focal, stand_size_pc2_drought = stand_size)

write.csv(pc2_drought_ses, "ses_PC2_drought.csv", row.names = FALSE)
