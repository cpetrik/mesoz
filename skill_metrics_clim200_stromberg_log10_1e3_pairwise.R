# Calculate different skill metrics for each ESM
# log10 transformed biomass (+ 1mgC/m2)
# compared against Stromberg model

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
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
ddir <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"
Amod <- read.csv(paste0(ddir,"climatol_All_hist_clims_200_stromberg_log10_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
Dmod <- read.csv(paste0(ddir,"climatol_DJF_hist_clims_200_stromberg_log10_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
Mmod <- read.csv(paste0(ddir,"climatol_MAM_hist_clims_200_stromberg_log10_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
Jmod <- read.csv(paste0(ddir,"climatol_JJA_hist_clims_200_stromberg_log10_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
Smod <- read.csv(paste0(ddir,"climatol_SON_hist_clims_200_stromberg_log10_1e4.csv"),sep=",",header = T,stringsAsFactors = F)

frval <- read.csv(paste0(ddir,"corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",header = T,stringsAsFactors = F)
frval <- frval[2:7,]
Kcor <- frval
Scor <- frval
Pcor <- frval

## Data exploration -------------------------------------------------------------------

## Correlations
for (x in 1:6) {
  amod <- as.data.frame(Amod$obsSM)
  amod[,2] <- Amod[,x+1]
  #amod <- na.omit(amod)
  Kcor[x,2] <- cor(amod[,1], amod[,2], method = "kendall", use="pairwise.complete.obs")
  
  dmod <- as.data.frame(Dmod$obsSM)
  dmod[,2] <- Dmod[,x+1]
  #dmod <- na.omit(dmod)
  Kcor[x,3] <- cor(dmod[,1], dmod[,2], method = "kendall", use="pairwise.complete.obs")
  
  mmod <- as.data.frame(Mmod$obsSM)
  mmod[,2] <- Mmod[,x+1]
  #mmod <- na.omit(mmod)
  Kcor[x,4] <- cor(mmod[,1], mmod[,2], method = "kendall", use="pairwise.complete.obs")
  
  jmod <- as.data.frame(Jmod$obsSM)
  jmod[,2] <- Jmod[,x+1]
  #jmod <- na.omit(jmod)
  Kcor[x,5] <- cor(jmod[,1], jmod[,2], method = "kendall", use="pairwise.complete.obs")
  
  smod <- as.data.frame(Smod$obsSM)
  smod[,2] <- Smod[,x+1]
  #smod <- na.omit(smod)
  Kcor[x,6] <- cor(smod[,1], smod[,2], method = "kendall", use="pairwise.complete.obs")
}

for (x in 1:6) {
  amod <- as.data.frame(Amod$obsSM)
  amod[,2] <- Amod[,x+1]
  #amod <- na.omit(amod)
  Scor[x,2] <- cor(amod[,1], amod[,2], method = "spearman", use="pairwise.complete.obs")
  
  dmod <- as.data.frame(Dmod$obsSM)
  dmod[,2] <- Dmod[,x+1]
  #dmod <- na.omit(dmod)
  Scor[x,3] <- cor(dmod[,1], dmod[,2], method = "spearman", use="pairwise.complete.obs")
  
  mmod <- as.data.frame(Mmod$obsSM)
  mmod[,2] <- Mmod[,x+1]
  #mmod <- na.omit(mmod)
  Scor[x,4] <- cor(mmod[,1], mmod[,2], method = "spearman", use="pairwise.complete.obs")
  
  jmod <- as.data.frame(Jmod$obsSM)
  jmod[,2] <- Jmod[,x+1]
  #jmod <- na.omit(jmod)
  Scor[x,5] <- cor(jmod[,1], jmod[,2], method = "spearman", use="pairwise.complete.obs")
  
  smod <- as.data.frame(Smod$obsSM)
  smod[,2] <- Smod[,x+1]
  #smod <- na.omit(smod)
  Scor[x,6] <- cor(smod[,1], smod[,2], method = "spearman", use="pairwise.complete.obs")
}

for (x in 1:6) {
  amod <- as.data.frame(Amod$obsSM)
  amod[,2] <- Amod[,x+1]
  #amod <- na.omit(amod)
  Pcor[x,2] <- cor(amod[,1], amod[,2], method = "pearson", use="pairwise.complete.obs")
  
  dmod <- as.data.frame(Dmod$obsSM)
  dmod[,2] <- Dmod[,x+1]
  #dmod <- na.omit(dmod)
  Pcor[x,3] <- cor(dmod[,1], dmod[,2], method = "pearson", use="pairwise.complete.obs")
  
  mmod <- as.data.frame(Mmod$obsSM)
  mmod[,2] <- Mmod[,x+1]
  #mmod <- na.omit(mmod)
  Pcor[x,4] <- cor(mmod[,1], mmod[,2], method = "pearson", use="pairwise.complete.obs")
  
  jmod <- as.data.frame(Jmod$obsSM)
  jmod[,2] <- Jmod[,x+1]
  #jmod <- na.omit(jmod)
  Pcor[x,5] <- cor(jmod[,1], jmod[,2], method = "pearson", use="pairwise.complete.obs")
  
  smod <- as.data.frame(Smod$obsSM)
  smod[,2] <- Smod[,x+1]
  #smod <- na.omit(smod)
  Pcor[x,6] <- cor(smod[,1], smod[,2], method = "pearson", use="pairwise.complete.obs")
}


write.table(Kcor,paste0(ddir,"Kendall_corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",row.names=T)
write.table(Scor,paste0(ddir,"Spearman_corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",row.names=T)
write.table(Pcor,paste0(ddir,"Pearson_corr_hist_clims_200_stromberg_log10_KK_1e4.csv"),sep=",",row.names=T)




