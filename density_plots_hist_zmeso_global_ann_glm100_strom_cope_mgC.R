# Density scatter plot of modeled zmeso biomass vs. obs
# By model-specific biomes
# GLMM100, COPEPOD, Stromberg models/obs

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
Zm <- read.csv(paste0(datap,"skill_hist_model_obsglm100_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
SZm <- read.csv(paste0(datap,"skill_hist_Stromberg_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
CZm <- read.csv(paste0(datap,"skill_hist_COPEPOD_all_clim_200_chls.csv"),sep=",",header = T,stringsAsFactors = F)

SzM <- as.data.frame(SZm[,c("Lat","Lon","Stromberg")])
CzM <- as.data.frame(CZm[,c("Lat","Lon","zmeso200")])
names(CzM) <- c("Lat","Lon","Copepod")

Zm1 <- Zm
Zm <- merge(Zm1,SzM,by=c("Lat","Lon"),all=TRUE)
Zm <- merge(Zm,CzM,by=c("Lat","Lon"),all=TRUE)

#Units from g to mg
Zm[,3:9] <- 1e3*Zm[,3:9]
Zm$Stromberg <- 1e3*Zm$Stromberg


### Plots -------------------------------------------------------
xlmts2 <- c( 1, 3e4 ) 

## Zoo ESMs vs Obs
c0 <- ggplot(Zm, aes(y=CAN, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("CAN") + xlab("") + ggtitle("obsGLMM") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m0 <- ggplot(Zm, aes(y=CMCC, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CMCC") + xlab("") + #ggtitle("") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n0 <- ggplot(Zm, aes(y=CNRM, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CNRM") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g0 <- ggplot(Zm, aes(y=GFDL, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 250) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("GFDL") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i0 <- ggplot(Zm, aes(y=IPSL, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("IPSL") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u0 <- ggplot(Zm, aes(y=UK, x=obsGLM)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("UK") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

cS <- ggplot(Zm, aes(y=CAN, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("obsSM") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
mS <- ggplot(Zm, aes(y=CMCC, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
nS <- ggplot(Zm, aes(y=CNRM, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
gS <- ggplot(Zm, aes(y=GFDL, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 250) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
iS <- ggplot(Zm, aes(y=IPSL, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
uS <- ggplot(Zm, aes(y=UK, x=Stromberg)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

cC <- ggplot(Zm, aes(y=CAN, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("COPEPOD") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
mC <- ggplot(Zm, aes(y=CMCC, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
nC <- ggplot(Zm, aes(y=CNRM, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
gC <- ggplot(Zm, aes(y=GFDL, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 250) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
iC <- ggplot(Zm, aes(y=IPSL, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
uC <- ggplot(Zm, aes(y=UK, x=Copepod)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')



## With Global
png(paste0(figp,'scatter_hist_zmeso_global_esm_vs_obs_log10_mgC.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 11*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c0,cS,cC,
           NULL,NULL,NULL,
           m0,mS,mC,
           NULL,NULL,NULL,
           n0,nS,nC,
           NULL,NULL,NULL,
           g0,gS,gC,
           NULL,NULL,NULL,
           i0,iS,iC,
           NULL,NULL,NULL,
           u0,uS,uC,
           nrow = 11, ncol = 3,
           rel_widths = c( 1, 1, 1 ), 
           rel_heights = c(1,-0.15,1,-0.15,1,-0.15,1,-0.15,1,-0.15,1) ,
           align = 'h' )
dev.off()



