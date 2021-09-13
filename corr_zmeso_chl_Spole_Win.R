# Calculate correlation of modeled zmeso biomass and chl

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
CM <- read.csv("can_hist_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_hist_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_hist_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_hist_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_hist_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Ann <- read.csv("obs_mod_chl_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## Winter DJF
aCM <- subset(CM, Lat < -60)
aNM <- subset(NM, Lat < -60)
aGM <- subset(GM, Lat < -60)
aIM <- subset(IM, Lat < -60)
aUM <- subset(UM, Lat < -60)
Amz <- subset(Ann, Lat < -60)

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chlDJF * 1e3
aGM$chl <- aGM$chlDJF * 1e3
aUM$chl <- aUM$chlDJF * 1e3

aNM$chl <- aNM$chlDJF 
aIM$chl <- aIM$chlDJF 

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,c("zooDJF","chl","sstDJF")])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,c("zooDJF","chl","sstDJF")])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,c("zooDJF","chl","sstDJF")])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,c("zooDJF","chl","sstDJF")])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,c("zooDJF","chl","sstDJF")])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)

Amydata <- as.matrix(Amz[,c("obs","chl")])
Amydata <- na.omit(Amydata)
Acorr <- rcorr(Amydata) #,type="pearson","spearman"
Acorr_all_p <- data.frame(Acorr$P)
Acorr_all_r <- data.frame(Acorr$r)

### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl
a1 <- ggplot(aCM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
   annotate( geom = 'text', y = 1e3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
   annotate( geom = 'text', y = 1e2, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e3, x = 2e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e3, x = 2e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 2e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 8e2, x = 2e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 2e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 8e2, x = 2e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 2e-5, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 2e-5, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

# untransform observations
Amz$Robs <- 10^Amz$obs
Amz$Rchl <- 10^Amz$chl * 1e-3 #from mg to g

a6 <- ggplot(Amz, aes(y=Robs, x=Rchl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("Obs zmeso mgC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 2e-5, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 2e-5, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_hist_mesoz_chl_Spole_Win_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,a6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sstDJF") + 
  annotate( geom = 'text', y = 75, x = -1, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 60, x = -1, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5) 

s2 <- ggplot(aNM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sstDJF") + 
  annotate( geom = 'text', y = 2000, x = -1, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1800, x = -1, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sstDJF") + 
  annotate( geom = 'text', y = 2800, x = -1, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2500, x = -1, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sstDJF") + 
  annotate( geom = 'text', y = 1200, x = -1, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1000, x = -1, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sstDJF") + 
  annotate( geom = 'text', y = 2000, x = -1, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1700, x = -1, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_hist_mesoz_sst_Spole_Win_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### LM of Chl -------------------------------------------------------------------
##Kwiatkowshki paper does not force through origin
# CAN
Clm <- lm(log10(zooDJF)~log10(chl), data=aCM)
summary(Clm)

Nlm <- lm(log10(zooDJF)~log10(chl), data=aNM)
summary(Nlm)

Glm <- lm(log10(zooDJF)~log10(chl), data=aGM)
summary(Glm)

Ilm <- lm(log10(zooDJF)~log10(chl), data=aIM)
summary(Ilm)

Ulm <- lm(log10(zooDJF)~log10(chl), data=aUM)
summary(Ulm)

Olm <- lm(log10(Robs)~log10(Rchl), data=Amz)
summary(Olm)

#SST
ClmS <- lm(zooDJF~sstDJF, data=aCM)
summary(ClmS)

NlmS <- lm(zooDJF~sstDJF, data=aNM)
summary(NlmS)

GlmS <- lm(zooDJF~sstDJF, data=aGM)
summary(GlmS)

IlmS <- lm(zooDJF~sstDJF, data=aIM)
summary(IlmS)

UlmS <- lm(zooDJF~sstDJF, data=aUM)
summary(UlmS)

## Save coefficients
Cff <- as.data.frame(t(coefficients(Clm)))
Cff[2,] <- t(coefficients(Nlm))
Cff[3,] <- t(coefficients(Glm))
Cff[4,] <- t(coefficients(Ilm))
Cff[5,] <- t(coefficients(Ulm))
Cff[6,] <- t(coefficients(Olm))
row.names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff,"Hist_coeffs_mod_chl_Spole_Win_clim_200.csv",sep=",",row.names=T)

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"Hist_coeffs_mod_sst_Spole_Win_200.csv",sep=",",row.names=T)


# =========================================================================
#  SSP 585 
# =========================================================================

# load data
CM <- read.csv("can_ssp585_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_ssp585_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_ssp585_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_ssp585_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_ssp585_mod_zoo_chl_sst_Win_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## Winter DJF
aCM <- subset(CM, Lat < -60)
aNM <- subset(NM, Lat < -60)
aGM <- subset(GM, Lat < -60)
aIM <- subset(IM, Lat < -60)
aUM <- subset(UM, Lat < -60)

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chlDJF * 1e3
aGM$chl <- aGM$chlDJF * 1e3
aUM$chl <- aUM$chlDJF * 1e3

aNM$chl <- aNM$chlDJF 
aIM$chl <- aIM$chlDJF 

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,c("zooDJF","chl","sstDJF")])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,c("zooDJF","chl","sstDJF")])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,c("zooDJF","chl","sstDJF")])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,c("zooDJF","chl","sstDJF")])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,c("zooDJF","chl","sstDJF")])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)



### Plots

## Zoo & Chl
a1 <- ggplot(aCM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e2, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 8e2, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 700, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zooDJF, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_chl_Spole_Win_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  annotate( geom = 'text', y = 100, x = -1, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = -1, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5)

s2 <- ggplot(aNM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  annotate( geom = 'text', y = 800, x = -1, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 700, x = -1, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  annotate( geom = 'text', y = 2600, x = -1, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2400, x = -1, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  annotate( geom = 'text', y = 1600, x = -1, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1400, x = -1, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zooDJF, x=sstDJF)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  annotate( geom = 'text', y = 2000, x = -1, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1800, x = -1, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_sst_Spole_Win_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### LM of Chl ------------------------------------------------------------------
##Kwiatkowshki paper does not force through origin, but are anomalies

# CAN
Clm <- lm(log10(zooDJF)~log10(chl), data=aCM)
summary(Clm)

Nlm <- lm(log10(zooDJF)~log10(chl), data=aNM)
summary(Nlm)

Glm <- lm(log10(zooDJF)~log10(chl), data=aGM)
summary(Glm)

Ilm <- lm(log10(zooDJF)~log10(chl), data=aIM)
summary(Ilm)

Ulm <- lm(log10(zooDJF)~log10(chl), data=aUM)
summary(Ulm)

#SST
ClmS <- lm(zooDJF~sstDJF, data=aCM)
summary(ClmS)

NlmS <- lm(zooDJF~sstDJF, data=aNM)
summary(NlmS)

GlmS <- lm(zooDJF~sstDJF, data=aGM)
summary(GlmS)

IlmS <- lm(zooDJF~sstDJF, data=aIM)
summary(IlmS)

UlmS <- lm(zooDJF~sstDJF, data=aUM)
summary(UlmS)


## Save coefficients
Cff <- as.data.frame(t(coefficients(Clm)))
Cff[2,] <- t(coefficients(Nlm))
Cff[3,] <- t(coefficients(Glm))
Cff[4,] <- t(coefficients(Ilm))
Cff[5,] <- t(coefficients(Ulm))
row.names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Cff,"SSP585_coeffs_mod_chl_Spole_Win_clim_200.csv",sep=",",row.names=T)

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"SSP585_coeffs_mod_sst_Spole_Win_200.csv",sep=",",row.names=T)

