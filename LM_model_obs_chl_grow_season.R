# Calculate correlation of obs and ESM chl

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
Chl <- read.csv("chl_obs_mod_mean_1950_2014_all_lats.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
## Latitude between 30 N and 30 S: all seasons
# All seasons
Amz <- subset(Chl, Lat > -30 & Lat < 30)

Chl <- na.omit(Chl)
Amz <- na.omit(Amz)

### Data exploration -------------------------------------------------------------------
## Correlations
Smydata <- as.matrix(Chl[,3:8])
Scorr <- rcorr(Smydata) #,type="pearson","spearman"
Scorr_all_p <- data.frame(Scorr$P)
Scorr_all_r <- data.frame(Scorr$r)

Amydata <- as.matrix(Amz[,3:8])
Acorr <- rcorr(Amydata) #,type="pearson","spearman"
Acorr_all_p <- data.frame(Acorr$P)
Acorr_all_r <- data.frame(Acorr$r)

## Use log axes instead
Ann <- Amz[,3:8]
Ann <- (10^Ann)
Ann <- as.data.frame(Ann)

Smm <- Chl[,3:8]
Smm <- (10^Smm)
Smm <- as.data.frame(Smm)

### Plots
xlmts2 <- c( 1e-6, 1e-2 ) 

## Annual
a1 <- ggplot(Ann, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(Ann, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[3,1],digits = 2)), size=5)

a3 <- ggplot(Ann, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[4,1],digits = 2)), size=5)

a4 <- ggplot(Ann, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[5,1],digits = 2)), size=5)

a5 <- ggplot(Ann, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[6,1],digits = 2)), size=5)

pdf( file = 'chl_model_obs_scatter_corr_tropics_only.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

## All lats
ylmts2 <- c( 3e-3, 10 ) 

m1 <- ggplot(Smm, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[2,1],digits = 2)), size=5)

m2 <- ggplot(Smm, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[3,1],digits = 2)), size=5)

m3 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[4,1],digits = 2)), size=5)

m4 <- ggplot(Smm, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[5,1],digits = 2)), size=5)

m5 <- ggplot(Smm, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-6, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-6, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[6,1],digits = 2)), size=5)

pdf( file = 'chl_model_obs_scatter_corr_all_lats.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( m1,m2,m3,m4,m5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### Use tropics only for LM
## May need to force through origin
#linear regression 
cmdl <- lm(Amz$CAN ~ Amz$obs+0)
cpred <- predict.lmc(mdl)
summary(cmdl)
#Amz$obs 1.030264
#Multiple R-squared:  0.979,	Adjusted R-squared:  0.979

nmdl <- lm(Amz$CNRM ~ Amz$obs +0)
npred <- predict.lmc(mdl)
summary(nmdl) 
#Amz$obs 0.9496320
#Multiple R-squared:  0.9956,	Adjusted R-squared:  0.9956 

gmdl <- lm(Amz$GFDL ~ Amz$obs+0)
gpred <- predict.lmc(mdl)
summary(gmdl) 
#Amz$obs 0.9188706
#Multiple R-squared:  0.995,	Adjusted R-squared:  0.995

imdl <- lm(Amz$IPSL ~ Amz$obs+0)
ipred <- predict.lmc(mdl)
summary(imdl) 
#Amz$obs 0.9442303
#Multiple R-squared:  0.9966,	Adjusted R-squared:  0.9966

umdl <- lm(Amz$UK ~ Amz$obs+0)
upred <- predict.lmc(mdl)
summary(umdl) 
#Amz$obs  0.9368357
#Multiple R-squared:  0.9935,	Adjusted R-squared:  0.9935

Cff <- as.data.frame(cmdl$coefficients)
Cff[,2] <- nmdl$coefficients
Cff[,3] <- gmdl$coefficients
Cff[,4] <- imdl$coefficients
Cff[,5] <- umdl$coefficients
names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UK")

##
Hdf <- Amz %>%
  select(obs,CAN,CNRM,GFDL,IPSL,UK) %>%
  gather(key = "model", value = "value", -obs)

b1 <- ggplot(Hdf, aes(x = obs, y = value)) + 
  geom_point(aes(color = model)) + 
  geom_smooth(aes(color = model), method="lm", se=FALSE, size = 0.75) + 
  ggtitle("Hist") + 
  ylab("log10 model chl (g m-3)") + xlab("log10 obs chl (g m-3)") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))


pdf(file = 'chl_model_obs_scatterLM_corr_tropics_only.pdf', width = unit( 10, 'cm' ), height = unit( 10, 'cm' ))
b1
dev.off()

write.table(Cff,"coeffs_mod_obs_log10chl_tropics_only.csv",sep=",",row.names=T)


