# Calculate different skill metrics for each ESM
# log transformed biomass

rm(list=ls())

library(ggplot2)
library(gridExtra)
library(cowplot) #plot_grid
library(RColorBrewer)
library("tidyverse") 

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Tr1 <- read.csv("skill_model_obsglm_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Tr2 <- read.csv("skill_model_obsglm_DJF_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Tr3 <- read.csv("skill_model_obsglm_MAM_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Tr4 <- read.csv("skill_model_obsglm_JJA_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Tr5 <- read.csv("skill_model_obsglm_SON_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

#----------------------------- Diff Data Scaling ------------------------------------------------
### Standardization 
#1st take log
Tr1[,3:8] <- log(Tr1[,3:8] + 1e-16)
Tr2[,3:8] <- log(Tr2[,3:8] + 1e-16)
Tr3[,3:8] <- log(Tr3[,3:8] + 1e-16)
Tr4[,3:8] <- log(Tr4[,3:8] + 1e-16)
Tr5[,3:8] <- log(Tr5[,3:8] + 1e-16)

#raw
Tl1 <- Tr1[,3:8]
Tl2 <- Tr2[,3:8]
Tl3 <- Tr3[,3:8]
Tl4 <- Tr4[,3:8]
Tl5 <- Tr5[,3:8]

# scaled -1 to 1
Ts1 <- apply(Tl1, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Ts1 <- as.data.frame(Ts1)
Ts2 <- apply(Tl2, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Ts2 <- as.data.frame(Ts2)
Ts3 <- apply(Tl3, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Ts3 <- as.data.frame(Ts3)
Ts4 <- apply(Tl4, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Ts4 <- as.data.frame(Ts4)
Ts5 <- apply(Tl5, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Ts5 <- as.data.frame(Ts5)

#----------------------------- Scatter plots ----------------------------- 
## Raw
#annual
a1 <- ggplot(Tl1, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CAN annual")
a2 <- ggplot(Tl1, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CNRM")
a3 <- ggplot(Tl1, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("GFDL")
a4 <- ggplot(Tl1, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("IPSL")
a5 <- ggplot(Tl1, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("UKESM")
#winter
d1 <- ggplot(Tl2, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CAN winter")
d2 <- ggplot(Tl2, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CNRM")
d3 <- ggplot(Tl2, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("GFDL")
d4 <- ggplot(Tl2, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("IPSL")
d5 <- ggplot(Tl2, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("UKESM")
#spring
m1 <- ggplot(Tl3, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CAN spring")
m2 <- ggplot(Tl3, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CNRM")
m3 <- ggplot(Tl3, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("GFDL")
m4 <- ggplot(Tl3, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("IPSL")
m5 <- ggplot(Tl3, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("UKESM")
#summer
j1 <- ggplot(Tl4, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CAN summer")
j2 <- ggplot(Tl4, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CNRM")
j3 <- ggplot(Tl4, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("GFDL")
j4 <- ggplot(Tl4, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("IPSL")
j5 <- ggplot(Tl4, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("UKESM")
#fall
s1 <- ggplot(Tl5, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CAN fall")
s2 <- ggplot(Tl4, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("CNRM")
s3 <- ggplot(Tl5, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("GFDL")
s4 <- ggplot(Tl5, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("IPSL")
s5 <- ggplot(Tl5, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("log mesozoo") + ggtitle("UKESM")


## Scaled
#annual
Sa1 <- ggplot(Ts1, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CAN annual")
Sa2 <- ggplot(Ts1, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CNRM")
Sa3 <- ggplot(Ts1, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("GFDL")
Sa4 <- ggplot(Ts1, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("IPSL")
Sa5 <- ggplot(Ts1, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("UKESM")
#winter
Sd1 <- ggplot(Ts2, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CAN winter")
Sd2 <- ggplot(Ts2, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CNRM")
Sd3 <- ggplot(Ts2, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("GFDL")
Sd4 <- ggplot(Ts2, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("IPSL")
Sd5 <- ggplot(Ts2, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("UKESM")
#spring
Sm1 <- ggplot(Ts3, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CAN spring")
Sm2 <- ggplot(Ts3, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CNRM")
Sm3 <- ggplot(Ts3, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("GFDL")
Sm4 <- ggplot(Ts3, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("IPSL")
Sm5 <- ggplot(Ts3, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("UKESM")
#summer
Sj1 <- ggplot(Ts4, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CAN summer")
Sj2 <- ggplot(Ts4, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CNRM")
Sj3 <- ggplot(Ts4, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("GFDL")
Sj4 <- ggplot(Ts4, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("IPSL")
Sj5 <- ggplot(Ts4, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("UKESM")
#fall
Ss1 <- ggplot(Ts5, aes(obsGLM, CAN)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CAN fall")
Ss2 <- ggplot(Ts4, aes(obsGLM, CNRM)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("CNRM")
Ss3 <- ggplot(Ts5, aes(obsGLM, GFDL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("GFDL")
Ss4 <- ggplot(Ts5, aes(obsGLM, IPSL)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("IPSL")
Ss5 <- ggplot(Ts5, aes(obsGLM, UK)) + geom_point(alpha=0.25) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") + 
  stat_smooth(method = "lm") + theme_bw() + 
  xlab("obs-GLMM") + ylab("scaled mesozoo") + ggtitle("UKESM")


## Multi-panel plots
# pdf( file = paste0(figp,'corr_scatter_hist_mesoz_log_obsglm.pdf'), width = unit( 7.5, 'in' ), height = unit( 10, 'in' ) )
# plot_grid( a1,a2,a3,a4,a5,
#            d1,d2,d3,d4,d5,
#            m1,m2,m3,m4,m5,
#            j1,j2,j3,j4,j5,
#            s1,s2,s3,s4,s5,
#            nrow = 5, ncol = 5,
#            rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
#            align = 'h' )
# dev.off()

png(paste0(figp,'corr_scatter_hist_mesoz_log_obsglm_alpha.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( a1,d1,m1,j1,s1,
           a2,d2,m2,j2,s2,
           a3,d3,m3,j3,s3,
           a4,d4,m4,j4,s4,
           a5,d5,m5,j5,s5,
           nrow = 5, ncol = 5,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

png(paste0(figp,'corr_scatter_hist_mesoz_neg11_obsglm_alpha.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( Sa1,Sd1,Sm1,Sj1,Ss1,
           Sa2,Sd2,Sm2,Sj2,Ss2,
           Sa3,Sd3,Sm3,Sj3,Ss3,
           Sa4,Sd4,Sm4,Sj4,Ss4,
           Sa5,Sd5,Sm5,Sj5,Ss5,
           nrow = 5, ncol = 5,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()
