# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
CM <- read.csv("can_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
### GET OBS CHL and SST w/BIOMES !
OM <- read.csv("obs_zmeso200_chl_sst_clim_biomes.csv",sep=",",header = T,stringsAsFactors = F)
OM <- na.omit(OM)

### Correct chl units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
CM$chl <- CM$chl * 1e3
GM$chl <- GM$chl * 1e3
UM$chl <- UM$chl * 1e3
OM$chl <- OM$chl * 1e-3

### Convert all zoo units to gC/m2
#all models in molC: 12.01 g C in 1 mol C
#obs in mgC
CM$zmeso <- CM$zmeso * 12.01
NM$zmeso <- NM$zmeso * 12.01
GM$zmeso <- GM$zmeso * 12.01
IM$zmeso <- IM$zmeso * 12.01
UM$zmeso <- UM$zmeso * 12.01
OM$zmeso <- OM$zmeso * 1e-3

CM$Lzmeso <- log(CM$zmeso+1e-16)
NM$Lzmeso <- log(NM$zmeso+1e-16)
GM$Lzmeso <- log(GM$zmeso+1e-16)
IM$Lzmeso <- log(IM$zmeso+1e-16)
UM$Lzmeso <- log(UM$zmeso+1e-16)
OM$Lzmeso <- log(OM$zmeso+1e-16)

CM$Lchl <- log(CM$chl+1e-16)
NM$Lchl <- log(NM$chl+1e-16)
GM$Lchl <- log(GM$chl+1e-16)
IM$Lchl <- log(IM$chl+1e-16)
UM$Lchl <- log(UM$chl+1e-16)
OM$Lchl <- log(OM$chl+1e-16)

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl

Ccorr <- ddply(CM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ncorr <- ddply(NM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Gcorr <- ddply(GM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Icorr <- ddply(IM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ucorr <- ddply(UM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ocorr <- ddply(OM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))


### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl all biomes together
a1 <- ggplot(CM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso gC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() #+ 
   # annotate( geom = 'text', y = 1e3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
   # annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(NM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("CNRM zmeso gC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a3 <- ggplot(GM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("GFDL zmeso gC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a4 <- ggplot(IM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("IPSL zmeso gC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a5 <- ggplot(UM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +
  ylab("UKESM zmeso gC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a6 <- ggplot(OM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("Obs zmeso gC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() 

pdf( file = paste0(figp,'corr_hist_mesoz_chl_biomes_all_log.pdf'), width = unit( 12, 'cm' ), height = unit( 10, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,a6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 1 LC
b1 <- ggplot(subset(CM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +  
  ylab("CAN zmeso gC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr[1,2],digits = 2)), size=5) 

b2 <- ggplot(subset(NM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("CNRM zmeso gC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr[1,2],digits = 2)), size=5) 

b3 <- ggplot(subset(GM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("GFDL zmeso gC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr[1,2],digits = 2)), size=5) 

b4 <- ggplot(subset(IM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("IPSL zmeso gC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr[1,2],digits = 2)), size=5) 

b5 <- ggplot(subset(UM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("UKESM zmeso gC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-6, hjust = 0, label=paste0("r = ",signif(Ucorr[1,2],digits = 2)), size=5) 

b6 <- ggplot(subset(OM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("Obs zmeso gC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ocorr[1,2],digits = 2)), size=5) 

pdf( file = paste0(figp,'corr_hist_mesoz_chl_biome1_LC_log.pdf'), width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( b1,b2,b3,b4,b5,b6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 2 HCSS
c1 <- ggplot(subset(CM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +  
  ylab("CAN zmeso gC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[2,2],digits = 2)), size=5) 

c2 <- ggplot(subset(NM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("CNRM zmeso gC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[2,2],digits = 2)), size=5) 

c3 <- ggplot(subset(GM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("GFDL zmeso gC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[2,2],digits = 2)), size=5) 

c4 <- ggplot(subset(IM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("IPSL zmeso gC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[2,2],digits = 2)), size=5) 

c5 <- ggplot(subset(UM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("UKESM zmeso gC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[2,2],digits = 2)), size=5) 

c6 <- ggplot(subset(OM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("Obs zmeso gC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ocorr[2,2],digits = 2)), size=5) 

pdf( file = paste0(figp,'corr_hist_mesoz_chl_biome2_HCSS_log.pdf'), width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( c1,c2,c3,c4,c5,c6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 3 HCPS
d1 <- ggplot(subset(CM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso gC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[3,2],digits = 2)), size=5) 

d2 <- ggplot(subset(NM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("CNRM zmeso gC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[3,2],digits = 2)), size=5) 

d3 <- ggplot(subset(GM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("GFDL zmeso gC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[3,2],digits = 2)), size=5) 

d4 <- ggplot(subset(IM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("IPSL zmeso gC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[3,2],digits = 2)), size=5) 

d5 <- ggplot(subset(UM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("UKESM zmeso gC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[3,2],digits = 2)), size=5) 

d6 <- ggplot(subset(OM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("Obs zmeso gC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ocorr[3,2],digits = 2)), size=5) 

pdf( file = paste0(figp,'corr_hist_mesoz_chl_biome3_HCPS_log.pdf'), width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( d1,d2,d3,d4,d5,d6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()



### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==1))
Clm2 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==2))
Clm3 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==3))
summary(Clm1)
summary(Clm2)
summary(Clm3)

#CNRM
Nlm1 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==1))
Nlm2 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==2))
Nlm3 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==3))
summary(Nlm1)
summary(Nlm2)
summary(Nlm3)

#GFDL
Glm1 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==1))
Glm2 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==2))
Glm3 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==3))
summary(Glm1)
summary(Glm2)
summary(Glm3)

#IPSL
Ilm1 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==1))
Ilm2 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==2))
Ilm3 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==3))
summary(Ilm1)
summary(Ilm2)
summary(Ilm3)

#UKESM
Ulm1 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==1))
Ulm2 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==2))
Ulm3 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==3))
summary(Ulm1)
summary(Ulm2)
summary(Ulm3)

#Obs
Olm1 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==1))
Olm2 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==2))
Olm3 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==3))
summary(Olm1)
summary(Olm2)
summary(Olm3)

## Save coefficients
Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Nlm1))
Cff1[3,] <- t(coefficients(Glm1))
Cff1[4,] <- t(coefficients(Ilm1))
Cff1[5,] <- t(coefficients(Ulm1))
Cff1[6,] <- t(coefficients(Olm1))
names(Cff1) <- c("Intercept","Lchl")
row.names(Cff1) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff1,"Hist_coeffs_mod_chl_biome1_LC.csv",sep=",",row.names=T)

Cff2 <- as.data.frame(t(coefficients(Clm2)))
Cff2[2,] <- t(coefficients(Nlm2))
Cff2[3,] <- t(coefficients(Glm2))
Cff2[4,] <- t(coefficients(Ilm2))
Cff2[5,] <- t(coefficients(Ulm2))
Cff2[6,] <- t(coefficients(Olm2))
names(Cff2) <- c("Intercept","Lchl")
row.names(Cff2) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff2,"Hist_coeffs_mod_chl_biome2_HCSS.csv",sep=",",row.names=T)

Cff3 <- as.data.frame(t(coefficients(Clm3)))
Cff3[2,] <- t(coefficients(Nlm3))
Cff3[3,] <- t(coefficients(Glm3))
Cff3[4,] <- t(coefficients(Ilm3))
Cff3[5,] <- t(coefficients(Ulm3))
Cff3[6,] <- t(coefficients(Olm3))
names(Cff3) <- c("Intercept","Lchl")
row.names(Cff3) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff3,"Hist_coeffs_mod_chl_biome3_HCPS.csv",sep=",",row.names=T)


### Save obs std error
pchl <- as.data.frame(logspace(-5,-2,20))
names(pchl) <- "Lchl"
Oerr1 <- as.data.frame(predict(Olm1, newdata = pchl, interval = "confidence"))
Oerr2 <- as.data.frame(predict(Olm2, newdata = pchl, interval = "confidence"))
Oerr3 <- as.data.frame(predict(Olm3, newdata = pchl, interval = "confidence"))
write.table(Oerr1,"Hist_stderr_predict_obs_chl_biome1.csv",sep=",",row.names=T)
write.table(Oerr2,"Hist_stderr_predict_obs_chl_biome2.csv",sep=",",row.names=T)
write.table(Oerr3,"Hist_stderr_predict_obs_chl_biome3.csv",sep=",",row.names=T)

Er <- as.data.frame(t(coef(summary(Olm1))[, 2]))
Er[2,] <- t(coef(summary(Olm2))[, 2])
Er[3,] <- t(coef(summary(Olm3))[, 2])
names(Er) <- c("Intercept","Lchl")
row.names(Er) <- c("1","2","3")
write.table(Er,"Hist_stderr_coeff_obs_chl_all_biomes.csv",sep=",",row.names=T)



