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
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Tr <- read.csv("skill_model_obsglm_SON_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

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

### --------------------------- Taylor diagrams -----------------------------
## V1
library(openair)

### Raw (log)
TrC <- Tr2[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- Tr2[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- Tr2[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- Tr2[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- Tr2[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_fall_200_raw.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t1 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4","blue","purple3","red3","maroon4"))
dev.off()


## Z score
TrC <- TrZ[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- TrZ[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- TrZ[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- TrZ[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- TrZ[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_fall_200_norm.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t2 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4", "blue", "purple3", "red3", "maroon4"))
dev.off()


## Anom
TrC <- TrA[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- TrA[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- TrA[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- TrA[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- TrA[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_fall_200_anom.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t3 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4", "blue", "purple3", "red3", "maroon4"))
dev.off()


## 0 to 1
TrC <- Tr01[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- Tr01[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- Tr01[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- Tr01[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- Tr01[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_fall_200_01.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t4 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4", "blue", "purple3", "red3", "maroon4"))
dev.off()


## -1 to 1
TrC <- Tr11[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- Tr11[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- Tr11[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- Tr11[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- Tr11[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_fall_200_neg11.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t5 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4", "blue", "purple3", "red3", "maroon4"))
dev.off()


# pdf( file = paste0(figp,'Taylor_obsglm_fall_200_all.pdf'), width = unit( 8, 'cm' ), height = unit( 10, 'cm' ) )
# plot_grid( t1,t2,t3,t4,t5,
#            nrow = 3, ncol = 2,
#            rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
#            align = 'h' )
# dev.off()
