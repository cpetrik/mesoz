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
Tr <- read.csv("obs_mod_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

Tr2 = na.omit(Tr) #3072 vs 3617


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
Sdev[2,1] <- (sd(Tr2$CAN)/sd(Tr2$obs))
Sdev[3,1] <- (sd(Tr2$CNRM)/sd(Tr2$obs))
Sdev[4,1] <- (sd(Tr2$GFDL)/sd(Tr2$obs))
Sdev[5,1] <- (sd(Tr2$IPSL)/sd(Tr2$obs))
Sdev[6,1] <- (sd(Tr2$UK)/sd(Tr2$obs))
#names(Sdev) <- names(Tr2)[3:8]
names(Sdev) <- "StdDev"

## RMSDs
RMSD <- as.data.frame(rmse(Tr2$obs,Tr2$obs))
RMSD[2,1] <- rmse(Tr2$obs,Tr2$CAN)
RMSD[3,1] <- rmse(Tr2$CNRM,Tr2$obs)
RMSD[4,1] <- rmse(Tr2$GFDL,Tr2$obs)
RMSD[5,1] <- rmse(Tr2$IPSL,Tr2$obs)
RMSD[6,1] <- rmse(Tr2$UK,Tr2$obs)
names(RMSD) <- "RMSD"


### ----------------------- Taylor diagrams --------------------------
library(openair)

TaylorDiagram(Tr2, obs = "obs", mod = "CAN")

TaylorDiagram(Tr2, obs = "obs", mod = "CNRM")

#,"GFDL","IPSL","UK"), group = "model")

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

Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

TaylorDiagram(Tr3, obs = "obs", mod = "esm", group = "model")

taylor.diagram(Tr2$obs, 
               Tr2$CAN, 
               add=FALSE,
               col="green",
               pch=1,
               pos.cor=TRUE,
               xlab="SD",
               ylab="RMSD",
               main="Taylor Diagram",
               show.gamma=TRUE,
               ngamma=3,
               sd.arcs=1,
               ref.sd=TRUE,
               grad.corr.lines=c(0.2,0.4,0.6,0.8,0.9),
               pcex=1,cex.axis=1,
               normalize=TRUE,
               mar=c(5,4,6,6),
               lwd=10,
               font=5,
               lty=3)


legend(6.5,7.5,cex=1.2,pt.cex=1.2,legend=c("CAN","CNRM","GFDL","IPSL","UK"),pch=c(1,2,3,4,5),col=c("green","blue","purple","red","brown"))
#legend(1.5,1.5,cex=1.2,pt.cex=1.2,legend=c("volcano"),pch=4,col=c("red"))


