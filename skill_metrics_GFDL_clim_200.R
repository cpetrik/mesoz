# Calculate different skill metrics for each ESM
# log transformed biomass

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
Ann <- read.csv("obs_mod_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Win <- read.csv("obs_mod_winter_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Spr <- read.csv("obs_mod_spring_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Smm <- read.csv("obs_mod_summer_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Fll <- read.csv("obs_mod_fall_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

Ann = na.omit(Ann) #10348 -> 9029
Win = na.omit(Win) #3617 -> 3072
Spr = na.omit(Spr) #4354 -> 3727
Smm = na.omit(Smm) #4729 -> 4022
Fll = na.omit(Fll) #3835 -> 3283

### Data exploration -------------------------------------------------------------------
## Correlations
Tmydata <- as.matrix(Ann[,3:8])
Rcorr <- rcorr(Tmydata) #,type="pearson","spearman"
Rcorr_all_p <- data.frame(Rcorr$P)
Rcorr_all_r <- data.frame(Rcorr$r)

Wmydata <- as.matrix(Win[,3:8])
Wcorr <- rcorr(Wmydata) #,type="pearson","spearman"
Wcorr_all_p <- data.frame(Wcorr$P)
Wcorr_all_r <- data.frame(Wcorr$r)

Smydata <- as.matrix(Spr[,3:8])
Scorr <- rcorr(Smydata) #,type="pearson","spearman"
Scorr_all_p <- data.frame(Scorr$P)
Scorr_all_r <- data.frame(Scorr$r)

Mmydata <- as.matrix(Smm[,3:8])
Mcorr <- rcorr(Mmydata) #,type="pearson","spearman"
Mcorr_all_p <- data.frame(Mcorr$P)
Mcorr_all_r <- data.frame(Mcorr$r)

Fmydata <- as.matrix(Fll[,3:8])
Fcorr <- rcorr(Fmydata) #,type="pearson","spearman"
Fcorr_all_p <- data.frame(Fcorr$P)
Fcorr_all_r <- data.frame(Fcorr$r)

## RMSDs
aRMSD <- as.data.frame(rmse(Ann$obs,Ann$obs))
aRMSD[2,1] <- rmse(Ann$obs,Ann$CAN)
aRMSD[3,1] <- rmse(Ann$CNRM,Ann$obs)
aRMSD[4,1] <- rmse(Ann$GFDL,Ann$obs)
aRMSD[5,1] <- rmse(Ann$IPSL,Ann$obs)
aRMSD[6,1] <- rmse(Ann$UK,Ann$obs)
names(aRMSD) <- "RMSD"

wRMSD <- as.data.frame(rmse(Win$obs,Win$obs))
wRMSD[2,1] <- rmse(Win$obs,Win$CAN)
wRMSD[3,1] <- rmse(Win$CNRM,Win$obs)
wRMSD[4,1] <- rmse(Win$GFDL,Win$obs)
wRMSD[5,1] <- rmse(Win$IPSL,Win$obs)
wRMSD[6,1] <- rmse(Win$UK,Win$obs)
names(wRMSD) <- "RMSD"

sRMSD <- as.data.frame(rmse(Spr$obs,Spr$obs))
sRMSD[2,1] <- rmse(Spr$obs,Spr$CAN)
sRMSD[3,1] <- rmse(Spr$CNRM,Spr$obs)
sRMSD[4,1] <- rmse(Spr$GFDL,Spr$obs)
sRMSD[5,1] <- rmse(Spr$IPSL,Spr$obs)
sRMSD[6,1] <- rmse(Spr$UK,Spr$obs)
names(sRMSD) <- "RMSD"

mRMSD <- as.data.frame(rmse(Smm$obs,Smm$obs))
mRMSD[2,1] <- rmse(Smm$obs,Smm$CAN)
mRMSD[3,1] <- rmse(Smm$CNRM,Smm$obs)
mRMSD[4,1] <- rmse(Smm$GFDL,Smm$obs)
mRMSD[5,1] <- rmse(Smm$IPSL,Smm$obs)
mRMSD[6,1] <- rmse(Smm$UK,Smm$obs)
names(mRMSD) <- "RMSD"

fRMSD <- as.data.frame(rmse(Fll$obs,Fll$obs))
fRMSD[2,1] <- rmse(Fll$obs,Fll$CAN)
fRMSD[3,1] <- rmse(Fll$CNRM,Fll$obs)
fRMSD[4,1] <- rmse(Fll$GFDL,Fll$obs)
fRMSD[5,1] <- rmse(Fll$IPSL,Fll$obs)
fRMSD[6,1] <- rmse(Fll$UK,Fll$obs)
names(fRMSD) <- "RMSD"

## Use log axes instead
Ann <- Ann[,3:8]
Ann <- (10^Ann)
Ann <- as.data.frame(Ann)

Win <- Win[,3:8]
Win <- (10^Win)
Win <- as.data.frame(Win)

Spr <- Spr[,3:8]
Spr <- (10^Spr)
Spr <- as.data.frame(Spr)

Smm <- Smm[,3:8]
Smm <- (10^Smm)
Smm <- as.data.frame(Smm)

Fll <- Fll[,3:8]
Fll <- (10^Fll)
Fll <- as.data.frame(Fll)

### Plots

## Annual
xlmts2 <- c( 1, 3e4 ) 

g1 <- ggplot(Ann, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[2,1],digits = 2)), size=5)

g2 <- ggplot(Ann, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

g3 <- ggplot(Ann, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[4,1],digits = 2)), size=5)

g4 <- ggplot(Ann, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[5,1],digits = 2)), size=5)

g5 <- ggplot(Ann, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.1, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_all_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( g1,g2,g3,g4,g5,
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
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[3,1],digits = 2)), size=5)

w3 <- ggplot(Win, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[4,1],digits = 2)), size=5)

w4 <- ggplot(Win, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[5,1],digits = 2)), size=5)

w5 <- ggplot(Win, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_winter_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
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
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[2,1],digits = 2)), size=5)

s2 <- ggplot(Spr, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[3,1],digits = 2)), size=5)

s3 <- ggplot(Spr, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[4,1],digits = 2)), size=5)

s4 <- ggplot(Spr, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[5,1],digits = 2)), size=5)

s5 <- ggplot(Spr, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_spring_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
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
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[2,1],digits = 2)), size=5)

m2 <- ggplot(Smm, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[3,1],digits = 2)), size=5)

m3 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[4,1],digits = 2)), size=5)

m4 <- ggplot(Smm, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[5,1],digits = 2)), size=5)

m5 <- ggplot(Smm, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_summer_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
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
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[2,1],digits = 2)), size=5)

f2 <- ggplot(Fll, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1.2e2, x = 3e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[3,1],digits = 2)), size=5)

f3 <- ggplot(Fll, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[4,1],digits = 2)), size=5)

f4 <- ggplot(Fll, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 50, x = 3e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 3e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[5,1],digits = 2)), size=5)

f5 <- ggplot(Fll, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1, x = 2e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_fall_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( f1,f2,f3,f4,f5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off


### GFDL only
xlmts <- c( 1, 1e5 ) 

c1 <- ggplot(Ann, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("Annual") +
  scale_y_log10() + scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[4,1],digits = 2)), size=5)

c2 <- ggplot(Win, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2")  + ggtitle("DJF") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Wcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(wRMSD[4,1],digits = 2)), size=5)

c3 <- ggplot(Spr, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("MAM") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 350, x = 5e3, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 250, x = 5e3, hjust = 0, label=paste0("rmse = ",signif(sRMSD[4,1],digits = 2)), size=5)

c4 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("JJA") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Mcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(mRMSD[4,1],digits = 2)), size=5)

c5 <- ggplot(Fll, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + ggtitle("SON") + scale_y_log10() + 
  scale_x_log10(limits = xlmts) + 
  annotate( geom = 'text', y = 190, x = 4e3, hjust = 0, label=paste0("r = ",signif(Fcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 4e3, hjust = 0, label=paste0("rmse = ",signif(fRMSD[4,1],digits = 2)), size=5)

pdf( file = 'corr_GFDL_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( c2,c3,c4,c5,c1,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off
