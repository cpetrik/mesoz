# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes
# Stromberg data only

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
#Zm <- read.csv(paste0(datap,"skill_hist_model_obsglm100_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
#Chl <- read.csv(paste0(datap,"hist_chl_model_obsglm100_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)
Zm <- read.csv(paste0(datap,"skill_hist_Stromberg_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv(paste0(datap,"hist_chl_seawifs_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)

SM <- as.data.frame(Zm)
SM[,5] <- as.data.frame(Chl[,3])
names(SM) <- c("Lat","Lon","mesoz","biome","chl")
SM <- na.omit(SM)

# OG <- as.data.frame(Zm[,c(1:2,9,16)])
# OG[,5] <- as.data.frame(Chl[,9])
# OG[,6] <- as.data.frame(Sst[,9])
# names(OG) <- c("Lat","Lon","mesoz","biome","chl","sst")
# OG <- na.omit(OG)

SM$Lzmeso <- log10(SM$mesoz+1e-16)
#OG$Lzmeso <- log10(OG$mesoz+1e-16)

SM$Lchl <- log10(SM$chl+1e-16)
#OG$Lchl <- log10(OG$chl+1e-16)

## Restrict latitudes to those observed
# maxlat <- max(OG$Lat)
# minlat <- min(OG$Lat)
# SM <- subset(SM, Lat >= minlat & Lat <= maxlat)

## Remove LC >45 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
uid1 <- which(SM$biome==1)
uid2 <- which(SM$Lat > 45 | SM$Lat < -45)
uid <- intersect(uid1,uid2)
SM <- SM[-uid,]

# oid1 <- which(OG$biome==1)
# oid2 <- which(OG$Lat > 45 | OG$Lat < -45)
# oid <- intersect(oid1,oid2)
# OG <- OG[-oid,]


### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#Stromberg 
Slm0 <- lm(Lzmeso~Lchl, data=SM)
Slm1 <- lm(Lzmeso~Lchl, data=subset(SM, biome==1))
Slm2 <- lm(Lzmeso~Lchl, data=subset(SM, biome==2))
Slm3 <- lm(Lzmeso~Lchl, data=subset(SM, biome==3))
summary(Slm1)
summary(Slm2)
summary(Slm3)

#Obs
# Olm0 <- lm(Lzmeso~Lchl, data=OG)
# Olm1 <- lm(Lzmeso~Lchl, data=subset(OG, biome==1))
# Olm2 <- lm(Lzmeso~Lchl, data=subset(OG, biome==2))
# Olm3 <- lm(Lzmeso~Lchl, data=subset(OG, biome==3))
# summary(Olm1)
# summary(Olm2)
# summary(Olm3)

## Save coefficients

Cff <- as.data.frame(t(coefficients(Slm0)))
Cff[2,] <- t(coefficients(Slm1))
Cff[3,] <- t(coefficients(Slm2))
Cff[4,] <- t(coefficients(Slm3))
names(Cff) <- c("Intercept","Lchl")
row.names(Cff) <- c("global","LC","HCSS","HCPS")
write.table(Cff,"Hist_coeffs_mod_chl_global_biomes_stromberg.csv",sep=",",row.names=T)


## Save standard error of obsGLMM
se <- as.data.frame(t(coefficients(Slm0)))
se[2,] <- t(coefficients(Slm1))
se[3,] <- t(coefficients(Slm2))
se[4,] <- t(coefficients(Slm3))

se[1,c(3:4)] <- as.data.frame(t(sqrt(diag(vcov(Slm0)))))
se[2,c(3:4)] <- sqrt(diag(vcov(Slm1)))
se[3,c(3:4)] <- sqrt(diag(vcov(Slm2)))
se[4,c(3:4)] <- sqrt(diag(vcov(Slm3)))
names(se) <- c("aSM","bSM","aSMse","bSMse")
row.names(se) <- c("Global","LC","HCSS","HCPS")
write.table(se,"Hist_std_err_biomes_stromberg.csv",sep=",",row.names=T)

## prediction tests
cmin <- log10(0.01)
cmax <- log10(25)
pchl <- exp(log(10)*seq(cmin,cmax,length=50))

plm0 <- predict(Slm0, newdata = data.frame(Lchl = log10(pchl)),
        interval = "confidence")
plm1 <- predict(Slm1, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm2 <- predict(Slm2, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm3 <- predict(Slm3, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")

CI0 <- as.data.frame(plm0)
CI1 <- as.data.frame(plm1)
CI2 <- as.data.frame(plm2)
CI3 <- as.data.frame(plm3)
CI0$Lchl <- log10(pchl)
CI1$Lchl <- log10(pchl)
CI2$Lchl <- log10(pchl)
CI3$Lchl <- log10(pchl)

write.table(CI0,"Hist_CIs_chl_global_stromberg.csv",sep=",",row.names=F)
write.table(CI1,"Hist_CIs_chl_biome1_LC_stromberg.csv",sep=",",row.names=F)
write.table(CI2,"Hist_CIs_chl_biome2_HCSS_stromberg.csv",sep=",",row.names=F)
write.table(CI3,"Hist_CIs_chl_biome3_HCPS_stromberg.csv",sep=",",row.names=F)
