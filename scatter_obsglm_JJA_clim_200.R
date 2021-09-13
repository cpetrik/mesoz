# Calculate different skill metrics for each ESM
# log transformed biomass

rm(list=ls())

#library(sm)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(PerformanceAnalytics)
library(Hmisc) #rcorr
#library(plyr)
library(cowplot) #plot_grid
library(RColorBrewer)
library(hydroGOF) #rmse
library("tidyverse") 

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Tr <- read.csv("skill_model_obsglm_JJA_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

#----------------------------- Diff Data Scaling ------------------------------------------------
### Standardization 
#1st take log
Tr[,3:8] <- log(Tr[,3:8] + 1e-16)
Tr2 <- Tr[,3:8]

# z-score standardize
TrZ <- as.data.frame(scale(Tr2))

# anomoly
TrA <- apply(Tr2, MARGIN = 2, FUN = function(X) (X - mean(X,na.rm=T)))
TrA <- as.data.frame(TrA)

# 0 to 1
Tr01 <- apply(Tr2, MARGIN = 2, FUN = function(X) (X - min(X,na.rm=T))/diff(range(X,na.rm=T)))
Tr01 <- as.data.frame(Tr01)

# -1 to 1
Tr11 <- apply(Tr2, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Tr11 <- as.data.frame(Tr11)


#----------------------------- Scatter plot data ----------------------------- 
### Raw
## Correlations
Tmydata <- as.matrix(Tr2)
Rcorr <- rcorr(Tmydata) #,type="pearson","spearman"
Rcorr_all_p <- data.frame(Rcorr$P)
Rcorr_all_r <- data.frame(Rcorr$r)

## RMSDs
rRMSD <- as.data.frame(rmse(Tr2$obsGLM,Tr2$obsGLM))
rRMSD[2,1] <- rmse(Tr2$CAN,Tr2$obsGLM,na.rm=T)
rRMSD[3,1] <- rmse(Tr2$CNRM,Tr2$obsGLM,na.rm=T)
rRMSD[4,1] <- rmse(Tr2$GFDL,Tr2$obsGLM,na.rm=T)
rRMSD[5,1] <- rmse(Tr2$IPSL,Tr2$obsGLM,na.rm=T)
rRMSD[6,1] <- rmse(Tr2$UK,Tr2$obsGLM,na.rm=T)
names(rRMSD) <- "RMSD"

### Group together with model as factor for facet wrap
Rdf1 <- Tr2 %>%
  select(CAN,CNRM,GFDL,IPSL,UK,obsGLM) %>%
  gather(key = "model", value = "mesoz", -obsGLM)
Ndf1 <- TrZ %>%
  select(CAN,CNRM,GFDL,IPSL,UK,obsGLM) %>%
  gather(key = "model", value = "mesoz", -obsGLM)
Adf1 <- TrA %>%
  select(CAN,CNRM,GFDL,IPSL,UK,obsGLM) %>%
  gather(key = "model", value = "mesoz", -obsGLM)
Zdf1 <- Tr01 %>%
  select(CAN,CNRM,GFDL,IPSL,UK,obsGLM) %>%
  gather(key = "model", value = "mesoz", -obsGLM)
Odf1 <- Tr11 %>%
  select(CAN,CNRM,GFDL,IPSL,UK,obsGLM) %>%
  gather(key = "model", value = "mesoz", -obsGLM)

#----------------------------- Scatter plots ----------------------------- 

p1 <- ggplot(Rdf1, aes(obsGLM, mesoz)) + 
  geom_point() + facet_wrap(. ~ model, nrow=3, ncol=2) + 
  stat_smooth(method = "lm") +
  theme_bw() + xlab("Summer obs-GLMM") + ylab("raw mesozooplankton")

p2 <- ggplot(Ndf1, aes(obsGLM, mesoz)) + 
  geom_point() + facet_wrap(. ~ model, nrow=3, ncol=2) + 
  stat_smooth(method = "lm") +
  theme_bw() + xlab("Summer obs-GLMM") + ylab("norm mesozooplankton")

p3 <- ggplot(Adf1, aes(obsGLM, mesoz)) + 
  geom_point() + facet_wrap(. ~ model, nrow=3, ncol=2) + 
  stat_smooth(method = "lm") +
  theme_bw() + xlab("Summer obs-GLMM") + ylab("anom mesozooplankton")

p4 <- ggplot(Zdf1, aes(obsGLM, mesoz)) + 
  geom_point() + facet_wrap(. ~ model, nrow=3, ncol=2) + 
  stat_smooth(method = "lm") +
  theme_bw() + xlab("Summer obs-GLMM") + ylab("0-1 mesozooplankton")

p5 <- ggplot(Odf1, aes(obsGLM, mesoz)) + 
  geom_point() + facet_wrap(. ~ model, nrow=3, ncol=2) + 
  stat_smooth(method = "lm") +
  theme_bw() + xlab("Summer obs-GLMM") + ylab("-1-1 mesozooplankton")

pdf( file = paste0(figp,'Scatter_wrap_obsglm_JJA_clim_200_raw.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
p1
dev.off()
pdf( file = paste0(figp,'Scatter_wrap_obsglm_JJA_clim_200_norm.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
p2
dev.off()
pdf( file = paste0(figp,'Scatter_wrap_obsglm_JJA_clim_200_anom.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
p3
dev.off()
pdf( file = paste0(figp,'Scatter_wrap_obsglm_JJA_clim_200_01.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
p4
dev.off()
pdf( file = paste0(figp,'Scatter_wrap_obsglm_JJA_clim_200_neg11.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
p5
dev.off()

