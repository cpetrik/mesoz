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

# load data
Tr <- read.csv("obs_mod_chl_stemp_falwin_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

#Tr2 <- na.omit(Tr) #3072 vs 3617
Tr2 <- Tr

### ----------------------- Taylor diagrams --------------------------
## V1
library(openair)

TaylorDiagram(Tr2, obs = "obs", mod = "CNRM")


TrC <- Tr2[,c("obs","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"

TrN <- Tr2[,c("obs","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"

TrG <- Tr2[,c("obs","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"

TrI <- Tr2[,c("obs","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"

TrU <- Tr2[,c("obs","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"


# With CAN
Tr3 <- rbind(TrN,TrG,TrI,TrU,TrC)

pdf( file = 'Taylor_stemp_clim_200_log.pdf') #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
TaylorDiagram(Tr3, obs = "obs", mod = "esm", group = "model",
              cols = c("darkorange2","red3","purple3", "blue", "green4"))
dev.off()


# Without CAN
NoC <- rbind(TrN,TrG,TrI,TrU)

pdf( file = 'Taylor_noCAN_stemp_clim_200_log.pdf') #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
TaylorDiagram(NoC, obs = "obs", mod = "esm", group = "model",
              cols = c("red3","purple3", "blue", "green4"))
dev.off()


