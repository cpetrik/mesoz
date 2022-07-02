# Calculate different skill metrics for each ESM
# log10 transformed biomass
# Kelly Kearney stats

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
frval <- read.csv(paste0(ddir,"corr_hist_clims_200_obsglm100_log10_KK.csv"),sep=",",header = T,stringsAsFactors = F)
frms <- read.csv(paste0(ddir,"urmse_hist_clims_200_obsglm100_log10_KK.csv"),sep=",",header = T,stringsAsFactors = F)
fbia <- read.csv(paste0(ddir,"bias_hist_clims_200_obsglm100_log10_KK.csv"),sep=",",header = T,stringsAsFactors = F)
fmae <- read.csv(paste0(ddir,"mae_hist_clims_200_obsglm100_log10_KK.csv"),sep=",",header = T,stringsAsFactors = F)

frval <- frval[2:7,] 
frms <- frms[2:7,] 
fbia <- fbia[2:7,] 
fmae <- fmae[2:7,] 

names(frval)[1] <- "Model"
names(frms)[1] <- "Model"
names(fbia)[1] <- "Model"
names(fmae)[1] <- "Model"

frval3 <- melt(frval)
names(frval3) <- c("Model","Season","Corr")
frms3 <- melt(frms)
names(frms3) <- c("Model","Season","RMSE")
fbia3 <- melt(fbia)
names(fbia3) <- c("Model","Season","Bias")
fmae3 <- melt(fmae)
names(fmae3) <- c("Model","Season","MAE")

frval3$Model <- factor(frval3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
frms3$Model <- factor(frms3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
fbia3$Model <- factor(fbia3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))
fmae3$Model <- factor(fmae3$Model,levels = c("UK","IPSL","GFDL","CNRM","CMCC","CAN"))

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
                       limit = c(0.3,0.05)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(RMSE,digits = 2)), color = "black", size = 3) 


frb <- ggplot(data = fbia3, aes(y=Model, x=Season, fill=Bias)) + 
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "PRGn", name="Bias", limit = c(-0.2,0.2)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(Bias,digits = 2)), color = "black", size = 3) 

fae <- ggplot(data = fmae3, aes(y=Model, x=Season, fill=MAE)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse", name="MAE", 
                       limit = c(0.3,0.04)) +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Season, Model, label = signif(MAE,digits = 2)), color = "black", size = 3) 


### ---------------------------- Together -----------------------
library(cowplot) #plot_grid

png(paste0(figp,'Heatmaps_urmse_clims_200_obsglm100_log10_KK.png'), 
    width = 4*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( fr,frm,frb,
           nrow = 3, ncol = 1,
           rel_widths = c(1,1), rel_heights = c(1,1) ,
           align = 'h' , labels = "auto", label_size = 12, hjust = -4)
dev.off()

png(paste0(figp,'Heatmaps_mae_clims_200_obsglm100_log10_KK.png'), 
    width = 4*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( fr,fae,frb,
           nrow = 3, ncol = 1,
           rel_widths = c(1,1), rel_heights = c(1,1) ,
           align = 'h' , labels = "auto", label_size = 12, hjust = -4)
dev.off()




