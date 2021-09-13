# Calculate correlation of modeled zmeso biomass and chl

rm(list=ls())

# library(sm)
# library(ggplot2)
# library(gridExtra)
# library(corrgram)
# library(PerformanceAnalytics)
library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid
# library(RColorBrewer)
# library(Metrics)


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
CM <- read.csv("can_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## All seasons
aCM <- subset(CM, Lat > -30 & Lat < 30)
aNM <- subset(NM, Lat > -30 & Lat < 30)
aGM <- subset(GM, Lat > -30 & Lat < 30)
aIM <- subset(IM, Lat > -30 & Lat < 30)
aUM <- subset(UM, Lat > -30 & Lat < 30)

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chl * 1e3
aGM$chl <- aGM$chl * 1e3
aUM$chl <- aUM$chl * 1e3

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,3:5])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,3:5])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,3:5])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,3:5])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,3:5])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)



### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl
a1 <- ggplot(aCM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl kg m-3") + 
  # annotate( geom = 'text', y = 600, x = 0, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 550, x = 0, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)
   scale_y_log10() + scale_x_log10() + 
   annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
   annotate( geom = 'text', y = 1e2, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 8e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 700, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_chl_tropics_ann_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  annotate( geom = 'text', y = 600, x = 15, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 550, x = 15, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5)
# scale_y_log10() + scale_x_log10() + 
# annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
# annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

s2 <- ggplot(aNM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  annotate( geom = 'text', y = 6000, x = 20, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 20, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  annotate( geom = 'text', y = 3300, x = 17, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3000, x = 17, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  annotate( geom = 'text', y = 1000, x = 15, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 900, x = 15, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  annotate( geom = 'text', y = 6000, x = 18, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 18, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_sst_tropics_ann_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### ============================= LM of Chl =================================
##Kwiatkowshki paper does not force through origin, but are anomalies

# CAN
Clm <- lm(log10(zoo)~log10(chl), data=aCM)
summary(Clm)

Nlm <- lm(log10(zoo)~log10(chl), data=aNM)
summary(Nlm)

Glm <- lm(log10(zoo)~log10(chl), data=aGM)
summary(Glm)

Ilm <- lm(log10(zoo)~log10(chl), data=aIM)
summary(Ilm)

Ulm <- lm(log10(zoo)~log10(chl), data=aUM)
summary(Ulm)


## Save coefficients
Cff <- as.data.frame(t(coefficients(Clm)))
Cff[2,] <- t(coefficients(Nlm))
Cff[3,] <- t(coefficients(Glm))
Cff[4,] <- t(coefficients(Ilm))
Cff[5,] <- t(coefficients(Ulm))
row.names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")

write.table(Cff,"SSP585_coeffs_mod_chl_trop_all_clim_200.csv",sep=",",row.names=T)






