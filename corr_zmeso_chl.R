# Calculate correlation of obs zmeso biomass and chl

rm(list=ls())

library(sm)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(PerformanceAnalytics)
library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid
library(RColorBrewer)
library(Metrics)


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
Ann <- read.csv("obs_mod_chl_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Win <- read.csv("obs_mod_chl_winter_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Spr <- read.csv("obs_mod_chl_spring_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Smm <- read.csv("obs_mod_chl_summer_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Fll <- read.csv("obs_mod_chl_fall_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
#Chl <- read.csv("seawifs_chl_ocx_growingseason_mean.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## Summer = JJA
Smz <- subset(Smm, Lat >= 60)

## Winter = DJF
Wmz <- subset(Win, Lat <= -60)

## Spring & Summer = JJA + MAM
Rmz <- subset(Spr, Lat < 60 & Lat >= 30)
Mmz <- subset(Smm, Lat < 60 & Lat >= 30)
#Combine
SSmz2 <- rbind(Rmz,Mmz)
#Use plyr to take mean of both
SSmz3 <- ddply( SSmz2, .( Lat,Lon ), summarize, 
             obs = mean( obs, na.rm = TRUE ), 
             CAN = mean( CAN, na.rm = TRUE ), 
             CNRM = mean( CNRM, na.rm = TRUE ), 
             GFDL = mean( GFDL, na.rm = TRUE ),
             IPSL = mean( IPSL, na.rm = TRUE ), 
             UK = mean( UK, na.rm = TRUE ), 
             chl = mean( chl, na.rm = TRUE ))
SSmz <- SSmz3

## Fall & Winter = SON + DJF
Fmz <- subset(Fll, Lat > -60 & Lat <= -30)
Nmz <- subset(Win, Lat > -60 & Lat <= -30)
#Combine
FWmz2 <- rbind(Fmz,Nmz)
#Use plyr to take mean of both
FWmz <- ddply( FWmz2, .( Lat,Lon ), summarize, 
                obs = mean( obs, na.rm = TRUE ), 
                CAN = mean( CAN, na.rm = TRUE ), 
                CNRM = mean( CNRM, na.rm = TRUE ), 
                GFDL = mean( GFDL, na.rm = TRUE ),
                IPSL = mean( IPSL, na.rm = TRUE ), 
                UK = mean( UK, na.rm = TRUE ), 
                chl = mean( chl, na.rm = TRUE ))

## All seasons
Amz <- subset(Ann, Lat > -30 & Lat < 30)

### Save for Taylor diagrams
write.table(Amz,"obs_mod_chl_trop_all_clim_200.csv",sep=",",row.names=F)
write.table(Wmz,"obs_mod_chl_spol_winter_clim_200.csv",sep=",",row.names=F)
write.table(Smz,"obs_mod_chl_npol_summer_clim_200.csv",sep=",",row.names=F)
write.table(SSmz,"obs_mod_chl_ntempl_sprsum_clim_200.csv",sep=",",row.names=F)
write.table(FWmz,"obs_mod_chl_stemp_falwin_clim_200.csv",sep=",",row.names=F)

### Data exploration -------------------------------------------------------------------
## Correlations
Smydata <- as.matrix(Smz[,3:9])
Scorr <- rcorr(Smydata) #,type="pearson","spearman"
Scorr_all_p <- data.frame(Scorr$P)
Scorr_all_r <- data.frame(Scorr$r)

Wmydata <- as.matrix(Wmz[,3:9])
Wcorr <- rcorr(Wmydata) #,type="pearson","spearman"
Wcorr_all_p <- data.frame(Wcorr$P)
Wcorr_all_r <- data.frame(Wcorr$r)

# Rmydata <- as.matrix(Rmz[,3:9])
# Rcorr <- rcorr(Rmydata) #,type="pearson","spearman"
# Rcorr_all_p <- data.frame(Rcorr$P)
# Rcorr_all_r <- data.frame(Rcorr$r)
# 
# Mmydata <- as.matrix(Mmz[,3:9])
# Mcorr <- rcorr(Mmydata) #,type="pearson","spearman"
# Mcorr_all_p <- data.frame(Mcorr$P)
# Mcorr_all_r <- data.frame(Mcorr$r)
#
# Fmydata <- as.matrix(Fmz[,3:9])
# Fcorr <- rcorr(Fmydata) #,type="pearson","spearman"
# Fcorr_all_p <- data.frame(Fcorr$P)
# Fcorr_all_r <- data.frame(Fcorr$r)
# 
# Nmydata <- as.matrix(Nmz[,3:9])
# Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
# Ncorr_all_p <- data.frame(Ncorr$P)
# Ncorr_all_r <- data.frame(Ncorr$r)

Amydata <- as.matrix(Amz[,3:9])
Acorr <- rcorr(Amydata) #,type="pearson","spearman"
Acorr_all_p <- data.frame(Acorr$P)
Acorr_all_r <- data.frame(Acorr$r)

FWmydata <- as.matrix(FWmz[,3:9])
FWcorr <- rcorr(FWmydata) #,type="pearson","spearman"
FWcorr_all_p <- data.frame(FWcorr$P)
FWcorr_all_r <- data.frame(FWcorr$r)

SSmydata <- as.matrix(SSmz[,3:9])
SScorr <- rcorr(SSmydata) #,type="pearson","spearman"
SScorr_all_p <- data.frame(SScorr$P)
SScorr_all_r <- data.frame(SScorr$r)


## RMSDs
library(hydroGOF)
aRMSD <- as.data.frame(rmse(Amz$obs,Amz$obs))
aRMSD[2,1] <- rmse(Amz$obs,Amz$CAN,na.rm = TRUE)
aRMSD[3,1] <- rmse(Amz$CNRM,Amz$obs,na.rm = TRUE)
aRMSD[4,1] <- rmse(Amz$GFDL,Amz$obs,na.rm = TRUE)
aRMSD[5,1] <- rmse(Amz$IPSL,Amz$obs,na.rm = TRUE)
aRMSD[6,1] <- rmse(Amz$UK,Amz$obs,na.rm = TRUE)
names(aRMSD) <- "RMSD"

wRMSD <- as.data.frame(rmse(Wmz$obs,Wmz$obs))
wRMSD[2,1] <- rmse(Wmz$obs,Wmz$CAN,na.rm = TRUE)
wRMSD[3,1] <- rmse(Wmz$CNRM,Wmz$obs,na.rm = TRUE)
wRMSD[4,1] <- rmse(Wmz$GFDL,Wmz$obs,na.rm = TRUE)
wRMSD[5,1] <- rmse(Wmz$IPSL,Wmz$obs,na.rm = TRUE)
wRMSD[6,1] <- rmse(Wmz$UK,Wmz$obs,na.rm = TRUE)
names(wRMSD) <- "RMSD"

sRMSD <- as.data.frame(rmse(SSmz$obs,SSmz$obs))
sRMSD[2,1] <- rmse(SSmz$obs,SSmz$CAN,na.rm = TRUE)
sRMSD[3,1] <- rmse(SSmz$CNRM,SSmz$obs,na.rm = TRUE)
sRMSD[4,1] <- rmse(SSmz$GFDL,SSmz$obs,na.rm = TRUE)
sRMSD[5,1] <- rmse(SSmz$IPSL,SSmz$obs,na.rm = TRUE)
sRMSD[6,1] <- rmse(SSmz$UK,SSmz$obs,na.rm = TRUE)
names(sRMSD) <- "RMSD"

mRMSD <- as.data.frame(rmse(Smz$obs,Smz$obs))
mRMSD[2,1] <- rmse(Smz$obs,Smz$CAN,na.rm = TRUE)
mRMSD[3,1] <- rmse(Smz$CNRM,Smz$obs,na.rm = TRUE)
mRMSD[4,1] <- rmse(Smz$GFDL,Smz$obs,na.rm = TRUE)
mRMSD[5,1] <- rmse(Smz$IPSL,Smz$obs,na.rm = TRUE)
mRMSD[6,1] <- rmse(Smz$UK,Smz$obs,na.rm = TRUE)
names(mRMSD) <- "RMSD"

fRMSD <- as.data.frame(rmse(FWmz$obs,FWmz$obs))
fRMSD[2,1] <- rmse(FWmz$obs,FWmz$CAN,na.rm = TRUE)
fRMSD[3,1] <- rmse(FWmz$CNRM,FWmz$obs,na.rm = TRUE)
fRMSD[4,1] <- rmse(FWmz$GFDL,FWmz$obs,na.rm = TRUE)
fRMSD[5,1] <- rmse(FWmz$IPSL,FWmz$obs,na.rm = TRUE)
fRMSD[6,1] <- rmse(FWmz$UK,FWmz$obs,na.rm = TRUE)
names(fRMSD) <- "RMSD"


## Use log axes instead
Ann <- Amz[,3:9]
Ann <- (10^Ann)
Ann <- as.data.frame(Ann)

Win <- Wmz[,3:9]
Win <- (10^Win)
Win <- as.data.frame(Win)

Spr <- SSmz[,3:9]
Spr <- (10^Spr)
Spr <- as.data.frame(Spr)

Smm <- Smz[,3:9]
Smm <- (10^Smm)
Smm <- as.data.frame(Smm)

Fll <- FWmz[,3:9]
Fll <- (10^Fll)
Fll <- as.data.frame(Fll)

### Plots

## Observations with Chlorophyll
xlmts2 <- c( 1, 3e4 ) 

g1 <- ggplot(Ann, aes(y=chl, x=obs)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("Tropics chl") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 10, x = 1, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[1,7],digits = 2)), size=5) + 
  annotate( geom = 'text', y = 5, x = 1, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[1,7],digits = 2)), size=5) 

g2 <- ggplot(Spr, aes(y=chl, x=obs)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("N Temp chl") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 10, x = 1, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[1,7],digits = 2)), size=5) + 
  annotate( geom = 'text', y = 5, x = 1, hjust = 0, label=paste0("p = ",signif(SScorr_all_p[1,7],digits = 2)), size=5) 

g3 <- ggplot(Smm, aes(y=chl, x=obs)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("N Polar chl") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 10, x = 1, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[1,7],digits = 2)), size=5) + 
  annotate( geom = 'text', y = 5, x = 1, hjust = 0, label=paste0("p = ",signif(SScorr_all_p[1,7],digits = 2)), size=5) 

g4 <- ggplot(Fll, aes(y=chl, x=obs)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("S Temp chl") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 2, x = 1, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[1,7],digits = 2)), size=5) + 
  annotate( geom = 'text', y = 1, x = 1, hjust = 0, label=paste0("p = ",signif(FWcorr_all_p[1,7],digits = 2)), size=5) 

g5 <- ggplot(Win, aes(y=chl, x=obs)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("S Polar") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 2, x = 1, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[1,7],digits = 2)), size=5) + 
  annotate( geom = 'text', y = 1, x = 1, hjust = 0, label=paste0("p = ",signif(Wcorr_all_p[1,7],digits = 2)), size=5) 

pdf( file = 'corr_chl_obs_grow_season_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( g2,g3,g4,g5,g1,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Annual
a1 <- ggplot(Ann, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[2,1],digits = 2)), size=5)

a2 <- ggplot(Ann, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

a3 <- ggplot(Ann, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[4,1],digits = 2)), size=5)

a4 <- ggplot(Ann, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[5,1],digits = 2)), size=5)

a5 <- ggplot(Ann, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_tropics_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Winter
ylmts2 <- c( 3e-3, 10 ) 

w1 <- ggplot(Win, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[2,1],digits = 2)), size=5)

w2 <- ggplot(Win, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[3,1],digits = 2)), size=5)

w3 <- ggplot(Win, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 300, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[4,1],digits = 2)), size=5)

w4 <- ggplot(Win, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 100, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[5,1],digits = 2)), size=5)

w5 <- ggplot(Win, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 20, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 15, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_Spolar_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( w1,w2,w3,w4,w5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Spring
ylmts2 <- c( 3e-3, 10 ) 

s1 <- ggplot(Spr, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[2,1],digits = 2)), size=5)

s2 <- ggplot(Spr, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[3,1],digits = 2)), size=5)

s3 <- ggplot(Spr, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[4,1],digits = 2)), size=5)

s4 <- ggplot(Spr, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[5,1],digits = 2)), size=5)

s5 <- ggplot(Spr, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-1, x = 2e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-2, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_Ntemp_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Summer
ylmts2 <- c( 3e-3, 10 ) 

m1 <- ggplot(Smm, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-5, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-6, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[2,1],digits = 2)), size=5)

m2 <- ggplot(Smm, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[3,1],digits = 2)), size=5)

m3 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[4,1],digits = 2)), size=5)

m4 <- ggplot(Smm, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[5,1],digits = 2)), size=5)

m5 <- ggplot(Smm, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 10, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_Npolar_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( m1,m2,m3,m4,m5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Fall
f1 <- ggplot(Fll, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-7, x = 2e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[2,1],digits = 2)), size=5)

f2 <- ggplot(Fll, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 300, x = 1e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 1e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[3,1],digits = 2)), size=5)

f3 <- ggplot(Fll, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 300, x = 2e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[4,1],digits = 2)), size=5)

f4 <- ggplot(Fll, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 300, x = 3e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[5,1],digits = 2)), size=5)

f5 <- ggplot(Fll, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 80, x = 2e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 50, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_Stemp_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( f1,f2,f3,f4,f5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


### GFDL only
xlmts <- c( 1, 1e5 ) 

c1 <- ggplot(Ann, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("Tropics") +
  scale_y_log10() + scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 130, x = 4e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 100, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[4,1],digits = 2)), size=5)

c2 <- ggplot(Win, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2")  + ggtitle("S Polar") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[4,1],digits = 2)), size=5)

c3 <- ggplot(Spr, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("N Temp") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 200, x = 5e3, hjust = 0, label=paste0("r = ",signif(SScorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 5e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[4,1],digits = 2)), size=5)

c4 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("N Polar") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[4,1],digits = 2)), size=5)

c5 <- ggplot(Fll, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("S Temp") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 300, x = 4e3, hjust = 0, label=paste0("r = ",signif(FWcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[4,1],digits = 2)), size=5)

pdf( file = 'corr_GFDL_lats_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( c2,c3,c4,c5,c1,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()
