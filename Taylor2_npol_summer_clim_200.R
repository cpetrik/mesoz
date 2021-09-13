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
Tr <- read.csv("obs_mod_chl_npol_summer_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

#Tr2 <- na.omit(Tr) #3072 vs 3617
Tr2 <- Tr

## Data exploration -------------------------------------------------------------------

### Skill metrica for Taylor diagram
# Taylor R
#num=nansum((o-omean).*(p-pmean));
#skill(7,j) = num/(n*osig*psig);

# Taylor normalized std
#skill(8,j) = psig/osig;

# total root mean square difference
#q1=nansum((o-p).^2);
#skill(10,j) = sqrt(q1/n);

## Correlations
Tmydata <- as.matrix(Tr2[,3:8])
Rcorr <- rcorr(Tmydata) #,type="pearson","spearman"
Rcorr_all_p <- data.frame(Rcorr$P)
Rcorr_all_r <- data.frame(Rcorr$r)

## Standard devs
Sdev <- as.data.frame(sd(Tr2$obs)/sd(Tr2$obs))
Sdev[2,1] <- (sd(Tr2$CAN,na.rm=T)/sd(Tr2$obs,na.rm=T))
Sdev[3,1] <- (sd(Tr2$CNRM,na.rm=T)/sd(Tr2$obs,na.rm=T))
Sdev[4,1] <- (sd(Tr2$GFDL,na.rm=T)/sd(Tr2$obs,na.rm=T))
Sdev[5,1] <- (sd(Tr2$IPSL,na.rm=T)/sd(Tr2$obs,na.rm=T))
Sdev[6,1] <- (sd(Tr2$UK,na.rm=T)/sd(Tr2$obs,na.rm=T))
#names(Sdev) <- names(Tr2)[3:8]
names(Sdev) <- "StdDev"

## RMSDs
library(hydroGOF)
RMSD <- as.data.frame(rmse(Tr2$obs,Tr2$obs))
RMSD[2,1] <- rmse(Tr2$obs,Tr2$CAN,na.rm=T)
RMSD[3,1] <- rmse(Tr2$CNRM,Tr2$obs,na.rm=T)
RMSD[4,1] <- rmse(Tr2$GFDL,Tr2$obs,na.rm=T)
RMSD[5,1] <- rmse(Tr2$IPSL,Tr2$obs,na.rm=T)
RMSD[6,1] <- rmse(Tr2$UK,Tr2$obs,na.rm=T)
names(RMSD) <- "RMSD"


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
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = 'Taylor_npol_clim_200_log.pdf') #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
TaylorDiagram(Tr3, obs = "obs", mod = "esm", group = "model",
              cols = c("darkorange2","red3","purple3", "blue", "green4"))
dev.off()


# Without CAN
NoC <- rbind(TrN,TrG,TrI,TrU)

pdf( file = 'Taylor_noCAN_npol_clim_200_log.pdf') #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
TaylorDiagram(NoC, obs = "obs", mod = "esm", group = "model",
              cols = c("red3","purple3", "blue", "green4"))
dev.off()


