# Calculate different skill metrics for each ESM
# 4th root transformed biomass

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

### -------------------------------- Raw 4rt -------------------------------
# load data
ddir <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"
frval <- read.csv(paste0(ddir,"corr_hist_clims_200_obsglm100_raw_4rt.csv"),sep=",",header = T,stringsAsFactors = F)
frms <- read.csv(paste0(ddir,"urmse_hist_clims_200_obsglm100_raw_4rt.csv"),sep=",",header = T,stringsAsFactors = F)
fbia <- read.csv(paste0(ddir,"bias_hist_clims_200_obsglm100_raw_4rt.csv"),sep=",",header = T,stringsAsFactors = F)

names(frval)[1] <- "Model"
names(frms)[1] <- "Model"
names(fbia)[1] <- "Model"

frval3 <- melt(frval)
names(frval3) <- c("Model","Season","Corr")
frms3 <- melt(frms)
names(frms3) <- c("Model","Season","RMSE")
fbia3 <- melt(fbia)
names(fbia3) <- c("Model","Season","Bias")

frval3$Model <- factor(frval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
frms3$Model <- factor(frms3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
fbia3$Model <- factor(fbia3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))

frms3$RMSE <- abs(frms3$RMSE)

### Heatmaps
fr <- ggplot(data = frval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.77,0.77),  
                       name="Pearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 


frm <- ggplot(data = frms3, aes(y=Model, x=Season, fill=RMSE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="RMSE", 
                       limit = c(0.37,0.11)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(RMSE,digits = 2)), color = "black", size = 3) 


frb <- ggplot(data = fbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.25,0.25)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 


### -------------------------------- -1 to 1 using 4th rt -------------------------------
# load data
zrval <- read.csv(paste0(ddir,"corr_hist_clims_200_obsglm100_neg11_4rt.csv"),sep=",",header = T,stringsAsFactors = F)
zrms <- read.csv(paste0(ddir,"urmse_hist_clims_200_obsglm100_neg11_4rt.csv"),sep=",",header = T,stringsAsFactors = F)
zbia <- read.csv(paste0(ddir,"bias_hist_clims_200_obsglm100_neg11_4rt.csv"),sep=",",header = T,stringsAsFactors = F)

names(zrval)[1] <- "Model"
names(zrms)[1] <- "Model"
names(zbia)[1] <- "Model"

zrval3 <- melt(zrval)
names(zrval3) <- c("Model","Season","Corr")
zrms3 <- melt(zrms)
names(zrms3) <- c("Model","Season","RMSE")
zbia3 <- melt(zbia)
names(zbia3) <- c("Model","Season","Bias")

zrval3$Model <- factor(zrval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zrms3$Model <- factor(zrms3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zbia3$Model <- factor(zbia3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))

zrms3$RMSE <- abs(zrms3$RMSE)

### Heatmaps
zr <- ggplot(data = zrval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "RdBu", limit = c(-0.77,0.77),  
                       name="Pearson\nCorrelation") +
  theme_minimal()+ labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

zrm <- ggplot(data = zrms3, aes(y=Model, x=Season, fill=RMSE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="RMSE", 
                       limit = c(1.7,0.6)) +
  theme_minimal() + labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(RMSE,digits = 2)), color = "black", size = 3) 

zrb <- ggplot(data = zbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.35,0.35)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + labs(y="")+
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 


### ---------------------------- Together -----------------------
library(cowplot) #plot_grid

png(paste0(figp,'Heatmaps_urmse_clims_200_obsglm100_raw_neg11_4rt.png'), 
    width = 9*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( fr,zr,frm,zrm,frb,zrb,
           nrow = 3, ncol = 2,
           rel_widths = c(1,1), rel_heights = c(1,1) ,
           align = 'h' , labels = "AUTO", label_size = 12, hjust = -4)
dev.off()




