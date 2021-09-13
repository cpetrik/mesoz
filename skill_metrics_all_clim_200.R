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

## cvF
#ylmts2 <- c( 3e-3, 10 ) 

p1 <- ggplot(Tr2, aes(y=CAN, x=obs)) + theme_bw(base_size=18) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN log10 zmeso") + xlab("log10 obs mgC m-2") + 
  annotate( geom = 'text', y = -8.1, x = 3.25, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = -8.7, x = 3.25, hjust = 0, label=paste0("rmse = ",signif(RMSD[2,1],digits = 2)), size=5)

p2 <- ggplot(Tr2, aes(y=CNRM, x=obs)) + theme_bw(base_size=18) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM log10 zmeso") + xlab("log10 obs mgC m-2") + 
  annotate( geom = 'text', y = 1.9, x = 3.25, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1.8, x = 3.25, hjust = 0, label=paste0("rmse = ",signif(RMSD[3,1],digits = 2)), size=5)

p3 <- ggplot(Tr2, aes(y=GFDL, x=obs)) + theme_bw(base_size=18) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL log10 zmeso") + xlab("log10 obs mgC m-2") + 
  annotate( geom = 'text', y = 2.2, x = 3.25, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2.1, x = 3.25, hjust = 0, label=paste0("rmse = ",signif(RMSD[4,1],digits = 2)), size=5)

p4 <- ggplot(Tr2, aes(y=IPSL, x=obs)) + theme_bw(base_size=18) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL log10 zmeso") + xlab("log10 obs mgC m-2") + 
  annotate( geom = 'text', y = 1.6, x = 3.25, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1.45, x = 3.25, hjust = 0, label=paste0("rmse = ",signif(RMSD[5,1],digits = 2)), size=5)

p5 <- ggplot(Tr2, aes(y=UK, x=obs)) + theme_bw(base_size=18) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK log10 zmeso") + xlab("log10 obs mgC m-2") + 
  annotate( geom = 'text', y = 0.3, x = 3.25, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 0.1, x = 3.25, hjust = 0, label=paste0("rmse = ",signif(RMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_all_clim_200.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( p1,p2,p3,p4,p5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Use log axes instead
Tr3 <- Tr2[,3:8]
Tr3 <- (10^Tr3)
Tr3 <- as.data.frame(Tr3)

ylmts2 <- c( 3e-3, 10 ) 

g1 <- ggplot(Tr3, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10() + 
  annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(RMSD[2,1],digits = 2)), size=5)

g2 <- ggplot(Tr3, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10() + 
  annotate( geom = 'text', y = 1.2e2, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 80, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(RMSD[3,1],digits = 2)), size=5)

g3 <- ggplot(Tr3, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10() + 
  annotate( geom = 'text', y = 180, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 150, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(RMSD[4,1],digits = 2)), size=5)

g4 <- ggplot(Tr3, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10() + 
  annotate( geom = 'text', y = 50, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 35, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(RMSD[5,1],digits = 2)), size=5)

g5 <- ggplot(Tr3, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK zmeso") + xlab("obs mgC m-2") + scale_y_log10() + 
  scale_x_log10() + 
  annotate( geom = 'text', y = 3, x = 2e3, hjust = 0, label=paste0("r = ",signif(Rcorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1.5, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(RMSD[6,1],digits = 2)), size=5)

pdf( file = 'corr_all_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( g1,g2,g3,g4,g5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### ----------------------- Taylor diagrams --------------------------
library(plotrix)

taylor.diagram(Tr2$obs, 
               Tr2$CAN, 
               add=FALSE,
               col="green",
               pch=1,
               pos.cor=TRUE,
               main="All climatology",
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

taylor.diagram(Tr2$obs, 
               Tr2$CNRM,
               add=TRUE,
               col="blue",
               pch=2,
               pos.cor=TRUE,
               main="All climatology",
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

taylor.diagram(Tr2$obs, 
               Tr2$GFDL, 
               add=TRUE,
               col="purple",
               pch=3,
               pos.cor=TRUE,
               main="All climatology",
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

taylor.diagram(Tr2$obs, 
               Tr2$IPSL,
               add=TRUE,
               col="red",
               pch=4,
               pos.cor=TRUE,
               xlab="SD",
               ylab="RMSD",
               main="All climatology",
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

taylor.diagram(Tr2$obs, 
               Tr2$UK,
               add=TRUE,
               col="brown",
               pch=5,
               pos.cor=TRUE,
               xlab="SD",
               ylab="RMSD",
               main="All climatology",
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

legend(6.75,7.75,cex=1.2,pt.cex=1.2,legend=c("CAN","CNRM","GFDL","IPSL","UK"),pch=c(1,2,3,4,5),col=c("green","blue","purple","red","brown"))

text(2,6.25,"Centered RMS Difference", col="darkgrey")


