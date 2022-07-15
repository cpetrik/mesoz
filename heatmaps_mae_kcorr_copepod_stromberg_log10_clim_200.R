# Calculate different skill metrics for each ESM
# log10 transformed biomass

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

### -------------------------------- COPEPOD -------------------------------
# load data
ddir <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"
frval <- read.csv(paste0(ddir,"corr_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
kval <- read.csv(paste0(ddir,"Kendall_corr_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
sval <- read.csv(paste0(ddir,"Spearman_corr_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
frms <- read.csv(paste0(ddir,"urmse_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
fbia <- read.csv(paste0(ddir,"bias_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
fmae <- read.csv(paste0(ddir,"mae_hist_clims_200_copepod_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)

frval <- frval[2:7,] 
frms <- frms[2:7,] 
fbia <- fbia[2:7,] 
fmae <- fmae[2:7,] 

names(frval)[1] <- "Model"
names(kval)[1] <- "Model"
names(sval)[1] <- "Model"
names(frms)[1] <- "Model"
names(fbia)[1] <- "Model"
names(fmae)[1] <- "Model"

frval3 <- melt(frval)
names(frval3) <- c("Model","Season","Corr")
kval3 <- melt(kval)
names(kval3) <- c("Model","Season","Corr")
sval3 <- melt(sval)
names(sval3) <- c("Model","Season","Corr")
frms3 <- melt(frms)
names(frms3) <- c("Model","Season","RMSE")
fbia3 <- melt(fbia)
names(fbia3) <- c("Model","Season","Bias")
fmae3 <- melt(fmae)
names(fmae3) <- c("Model","Season","MAE")

frval3$Model <- factor(frval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
kval3$Model <- factor(kval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
sval3$Model <- factor(sval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
frms3$Model <- factor(frms3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
fbia3$Model <- factor(fbia3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
fmae3$Model <- factor(fmae3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))

frms3$RMSE <- abs(frms3$RMSE)

### Heatmaps
frp <- ggplot(data = frval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),  
                       name="Pearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

frk <- ggplot(data = kval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="Kendall\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

frs <- ggplot(data = sval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Spearman\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 


frb <- ggplot(data = fbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.35,0.35)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 

fae <- ggplot(data = fmae3, aes(y=Model, x=Season, fill=MAE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="MAE", 
                       limit = c(0.6,0.33)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(MAE,digits = 2)), color = "black", size = 3) 


### -------------------------------- Stromberg -------------------------------
# load data
zrval <- read.csv(paste0(ddir,"corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
zkval <- read.csv(paste0(ddir,"Kendall_corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
zsval <- read.csv(paste0(ddir,"Spearman_corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
zrms <- read.csv(paste0(ddir,"urmse_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
zbia <- read.csv(paste0(ddir,"bias_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
zmae <- read.csv(paste0(ddir,"mae_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)

zrval <- zrval[2:7,] 
zrms <- zrms[2:7,] 
zbia <- zbia[2:7,] 
zmae <- zmae[2:7,] 

names(zrval)[1] <- "Model"
names(zkval)[1] <- "Model"
names(zsval)[1] <- "Model"
names(zrms)[1] <- "Model"
names(zbia)[1] <- "Model"
names(zmae)[1] <- "Model"

zrval3 <- melt(zrval)
names(zrval3) <- c("Model","Season","Corr")
zkval3 <- melt(zkval)
names(zkval3) <- c("Model","Season","Corr")
zsval3 <- melt(zsval)
names(zsval3) <- c("Model","Season","Corr")
zrms3 <- melt(zrms)
names(zrms3) <- c("Model","Season","RMSE")
zbia3 <- melt(zbia)
names(zbia3) <- c("Model","Season","Bias")
zmae3 <- melt(zmae)
names(zmae3) <- c("Model","Season","MAE")

zrval3$Model <- factor(zrval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zkval3$Model <- factor(zkval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zsval3$Model <- factor(zsval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zrms3$Model <- factor(zrms3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zbia3$Model <- factor(zbia3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
zmae3$Model <- factor(zmae3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))

zrms3$RMSE <- abs(zrms3$RMSE)

### Heatmaps
zrp <- ggplot(data = zrval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),  
                       name="Pearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

zrk <- ggplot(data = zkval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="Kendall\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

zrs <- ggplot(data = zsval3, aes(y=Model, x=Season, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Spearman\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Corr,digits = 2)), color = "black", size = 3) 

zrb <- ggplot(data = zbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.4,0.4)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 

zae <- ggplot(data = zmae3, aes(y=Model, x=Season, fill=MAE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="MAE", 
                       limit = c(0.51,0.23)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(MAE,digits = 2)), color = "black", size = 3) 


### ---------------------------- Together -----------------------
library(cowplot) #plot_grid

png(paste0(figp,'Heatmaps_mae_clims_200_copepod_stromberg_log10.png'), 
    width = 9*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( frk,zrk,fae,zae,frb,zrb,
           nrow = 3, ncol = 2,
           rel_widths = c(1,1), rel_heights = c(1,1) ,
           align = 'h' , labels = "auto", label_size = 12, hjust = -4)
dev.off()




