# Calculate different skill metrics for each ESM
# log transformed biomass

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
library(scico)
library(reshape2)

### -------------------------------- Raw -------------------------------
# load data
rval <- read.csv("corr_clims_200_obsglm_raw.csv",sep=",",header = T,stringsAsFactors = F)
rms <- read.csv("rmse_clims_200_obsglm_raw.csv",sep=",",header = T,stringsAsFactors = F)
bia <- read.csv("bias_clims_200_obsglm_raw.csv",sep=",",header = T,stringsAsFactors = F)

rval <- rval[2:5,]
rms <- rms[2:5,]
bia <- bia[2:5,]

names(rval)[1] <- "Model"
names(rms)[1] <- "Model"
names(bia)[1] <- "Model"

rval3 <- melt(rval)
names(rval3) <- c("Model","Season","Corr")
rms3 <- melt(rms)
names(rms3) <- c("Model","Season","RMSE")
bia3 <- melt(bia)
names(bia3) <- c("Model","Season","Bias")

rval3$Model <- factor(rval3$Model,levels = c("UK","IPSL","GFDL","CNRM"))
rms3$Model <- factor(rms3$Model,levels = c("UK","IPSL","GFDL","CNRM"))
bia3$Model <- factor(bia3$Model,levels = c("UK","IPSL","GFDL","CNRM"))

### Heatmaps
gr <- ggplot(data = rval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),  
                       name="Pearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 


grm <- ggplot(data = rms3, aes(y=Model, x=Season, fill=RMSE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="RMSE", 
                       limit = c(2.5,0.5)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(RMSE,digits = 2)), color = "black", size = 3) 


grb <- ggplot(data = bia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-1.5,1.5)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 



### -------------------------------- Raw -------------------------------
# load data
nrval <- read.csv("corr_clims_200_obsglm_neg11.csv",sep=",",header = T,stringsAsFactors = F)
nrms <- read.csv("rmse_clims_200_obsglm_neg11.csv",sep=",",header = T,stringsAsFactors = F)
nbia <- read.csv("bias_clims_200_obsglm_neg11.csv",sep=",",header = T,stringsAsFactors = F)

nrval <- nrval[2:5,]
nrms <- nrms[2:5,]
nbia <- nbia[2:5,]

names(nrval)[1] <- "Model"
names(nrms)[1] <- "Model"
names(nbia)[1] <- "Model"

nrval3 <- melt(nrval)
names(nrval3) <- c("Model","Season","Corr")
nrms3 <- melt(nrms)
names(nrms3) <- c("Model","Season","RMSE")
nbia3 <- melt(nbia)
names(nbia3) <- c("Model","Season","Bias")

nrval3$Model <- factor(nrval3$Model,levels = c("UK","IPSL","GFDL","CNRM"))
nrms3$Model <- factor(nrms3$Model,levels = c("UK","IPSL","GFDL","CNRM"))
nbia3$Model <- factor(nbia3$Model,levels = c("UK","IPSL","GFDL","CNRM"))


### Heatmaps
sr <- ggplot(data = nrval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),  
                       name="Pearson\nCorrelation") +
  theme_minimal()+ labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

srm <- ggplot(data = nrms3, aes(y=Model, x=Season, fill=RMSE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="RMSE", 
                       limit = c(0.65,0.2)) +
  theme_minimal() + labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(RMSE,digits = 2)), color = "black", size = 3) 

srb <- ggplot(data = nbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.65,0.65)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + labs(y="")+
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 


### ---------------------------- Together -----------------------
library(cowplot) #plot_grid

png(paste0(figp,'Heatmaps_clims_200_obsglm_raw_neg11_cb_noCAN.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( gr,sr,grm,srm,grb,srb,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()




