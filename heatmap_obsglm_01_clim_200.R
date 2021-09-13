# Calculate different skill metrics for each ESM
# log transformed biomass
# scaled 0 to 1

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"


library(Hmisc)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(scatterplot3d)
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)
library(gplots)
library("corrplot")

# load data
rval <- read.csv("corr_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)
nstd <- read.csv("nstd_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)
rms <- read.csv("rmse_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)
urms <- read.csv("urmse_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)
trms <- read.csv("trmse_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)
bia <- read.csv("bias_clims_200_obsglm_01.csv",sep=",",header = T,stringsAsFactors = F)

names(rval)[1] <- "Model"
names(nstd)[1] <- "Model"
names(rms)[1] <- "Model"
names(urms)[1] <- "Model"
names(trms)[1] <- "Model"
names(bia)[1] <- "Model"

### Heatmaps
#ggplot tile
library(reshape2)
rval3 <- melt(rval)
names(rval3) <- c("Model","Season","Corr")
  
gr <- ggplot(data = rval3, aes(x=Model, y=Season, fill=Corr)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(Corr,digits = 2)), color = "black", size = 4) 

png(paste0(figp,"Heatmap_corr_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gr
dev.off()

## RMSE
rms3 <- melt(rms)
names(rms3) <- c("Model","Season","RMSE")

grm <- ggplot(data = rms3, aes(x=Model, y=Season, fill=RMSE)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "magma", 
                       name="RMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(RMSE,digits = 2)), color = "white", size = 4) 

png(paste0(figp,"Heatmap_rmse_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
grm
dev.off()


## Unbiased RMSE
urms3 <- melt(urms)
names(urms3) <- c("Model","Season","unbRMSE")

gru <- ggplot(data = urms3, aes(x=Model, y=Season, fill=unbRMSE)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "magma", 
                       name="unbiased\nRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(unbRMSE,digits = 2)), color = "white", size = 4) 

png(paste0(figp,"Heatmap_urmse_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gru
dev.off()


## Total RMSE
trms3 <- melt(trms)
names(trms3) <- c("Model","Season","totRMSE")

grt <- ggplot(data = trms3, aes(x=Model, y=Season, fill=totRMSE)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "magma", 
                       name="total\nRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(totRMSE,digits = 2)), color = "white", size = 4) 

png(paste0(figp,"Heatmap_trmse_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
grt
dev.off()


## Bias
bia3 <- melt(bia)
names(bia3) <- c("Model","Season","Bias")

grb <- ggplot(data = bia3, aes(x=Model, y=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "magma", 
                       name="Bias") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(Bias,digits = 2)), color = "white", size = 4) 

png(paste0(figp,"Heatmap_bias_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
grb
dev.off()


## nstd
nstd3 <- melt(nstd)
names(nstd3) <- c("Model","Season","nstd")

grs <- ggplot(data = nstd3, aes(x=Model, y=Season, fill=nstd)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "magma", 
                       name="normalized\ns.d.") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Model, Season, label = signif(nstd,digits = 2)), color = "white", size = 4) 

png(paste0(figp,"Heatmap_nstd_clims_200_obsglm_01.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
grs
dev.off()
