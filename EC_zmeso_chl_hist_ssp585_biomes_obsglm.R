# Linear regressions of zmeso biomass with surf chl 
# By biome, annual climatology
# Historic and SSP585

rm(list=ls())

library(plyr)
library(cowplot) #plot_grid
library(R.matlab)
library(rdetools)
library("tidyverse")

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
Hist0 <- read.csv(paste0(datap,"Hist_coeffs_mod_chl_global_obsglm.csv"),sep=",",header = T,stringsAsFactors = F)
Hist1 <- read.csv(paste0(datap,"Hist_coeffs_mod_chl_biome1_LC45_obsglm.csv"),sep=",",header = T,stringsAsFactors = F)
Hist2 <- read.csv(paste0(datap,"Hist_coeffs_mod_chl_biome2_HCSS_obsglm.csv"),sep=",",header = T,stringsAsFactors = F)
Hist3 <- read.csv(paste0(datap,"Hist_coeffs_mod_chl_biome3_HCPS_obsglm.csv"),sep=",",header = T,stringsAsFactors = F)

SSP0 <- read.csv(paste0(datap,"SSP585_coeffs_mod_chl_global.csv"),sep=",",header = T,stringsAsFactors = F)
SSP1 <- read.csv(paste0(datap,"SSP585_coeffs_mod_chl_biome1_LC45.csv"),sep=",",header = T,stringsAsFactors = F)
SSP2 <- read.csv(paste0(datap,"SSP585_coeffs_mod_chl_biome2_HCSS.csv"),sep=",",header = T,stringsAsFactors = F)
SSP3 <- read.csv(paste0(datap,"SSP585_coeffs_mod_chl_biome3_HCPS.csv"),sep=",",header = T,stringsAsFactors = F)

## predict mesozoo over chl values
log10(0.01)
log10(25)
chl <- logspace(-2,2,20)
hzoo0 <- matrix(NA, nrow = 7, ncol = 20) #7 with obs
hzoo1 <- matrix(NA, nrow = 7, ncol = 20)
hzoo2 <- matrix(NA, nrow = 7, ncol = 20)
hzoo3 <- matrix(NA, nrow = 7, ncol = 20)
for(i in 1:7)
{
  hzoo0[i,] <- Hist0$Intercept[i] + Hist0$Lchl[i] * log10(chl)
  hzoo1[i,] <- Hist1$Intercept[i] + Hist1$Lchl[i] * log10(chl)
  hzoo2[i,] <- Hist2$Intercept[i] + Hist2$Lchl[i] * log10(chl)
  hzoo3[i,] <- Hist3$Intercept[i] + Hist3$Lchl[i] * log10(chl)
}

szoo0 <- matrix(NA, nrow = 6, ncol = 20)
szoo1 <- matrix(NA, nrow = 6, ncol = 20)
szoo2 <- matrix(NA, nrow = 6, ncol = 20)
szoo3 <- matrix(NA, nrow = 6, ncol = 20)
for(i in 1:6)
{
  szoo0[i,] <- SSP0$Intercept[i] + SSP0$Lchl[i] * log10(chl)
  szoo1[i,] <- SSP1$Intercept[i] + SSP1$Lchl[i] * log10(chl)
  szoo2[i,] <- SSP2$Intercept[i] + SSP2$Lchl[i] * log10(chl)
  szoo3[i,] <- SSP3$Intercept[i] + SSP3$Lchl[i] * log10(chl)
}

hzoo0 <- as.data.frame(t(hzoo0))
hzoo1 <- as.data.frame(t(hzoo1))
hzoo2 <- as.data.frame(t(hzoo2))
hzoo3 <- as.data.frame(t(hzoo3))
names(hzoo0) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo1) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo2) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo3) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")

szoo0 <- as.data.frame(t(szoo0))
szoo1 <- as.data.frame(t(szoo1))
szoo2 <- as.data.frame(t(szoo2))
szoo3 <- as.data.frame(t(szoo3))
names(szoo0) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")
names(szoo1) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")
names(szoo2) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")
names(szoo3) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

hzoo0$chl <- chl
hzoo1$chl <- chl
hzoo2$chl <- chl
hzoo3$chl <- chl
szoo0$chl <- chl
szoo1$chl <- chl
szoo2$chl <- chl
szoo3$chl <- chl


Hdf0 <- hzoo0 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM,obs) %>%
  gather(key = "model", value = "value", -chl)
Hdf1 <- hzoo1 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM,obs) %>%
  gather(key = "model", value = "value", -chl)
Hdf2 <- hzoo2 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM,obs) %>%
  gather(key = "model", value = "value", -chl)
Hdf3 <- hzoo3 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM,obs) %>%
  gather(key = "model", value = "value", -chl)

Sdf0 <- szoo0 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)
Sdf1 <- szoo1 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)
Sdf2 <- szoo2 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)
Sdf3 <- szoo3 %>%
  select(chl,CAN,CMCC,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)

## figures
ylmts0 <- c(-10,5)
h0 <- ggplot(Hdf0, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("Global Hist") + ylim(ylmts0) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

s0 <- ggplot(Sdf0, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("Global SSP 585") + ylim(ylmts0) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))


ylmts <- c(-10,10)
h1 <- ggplot(Hdf1, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + ylim(ylmts) +
  scale_x_log10() + ggtitle("LC Hist") + 
  ylab("log10 zmeso mgC m-2") + xlab("log10 chl mg m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

s1 <- ggplot(Sdf1, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + ylim(ylmts) +
  scale_x_log10() + ggtitle("LC SSP 585") + 
  ylab("log10 zmeso mgC m-2") + xlab("log10 chl mg m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

ylmts2 <- c(-5,5)
h2 <- ggplot(Hdf2, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCSS Hist") + ylim(ylmts2) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

s2 <- ggplot(Sdf2, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCSS SSP 585") + ylim(ylmts2) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

ylmts3 <- c(-5,5)
h3 <- ggplot(Hdf3, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCPS Hist") + ylim(ylmts3) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))

s3 <- ggplot(Sdf3, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCPS SSP 585") + ylim(ylmts3) +
  ylab("log10 zmeso gC m-2") + xlab("log10 chl g m-3") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#555555","#AA3377"))


### linear regression of hist vs. future -----------------------------
mdl0 <- lm(SSP0$Lchl ~ Hist0$Lchl[1:6])
ypred0 <- predict.lm(mdl0)
summary(mdl0)

mdl0_noC <- lm(SSP0$Lchl[c(2:6)] ~ Hist0$Lchl[c(2:6)])
ypred0_noC <- predict.lm(mdl0_noC)
summary(mdl0_noC)

mdl1 <- lm(SSP1$Lchl ~ Hist1$Lchl[1:6])
ypred1 <- predict.lm(mdl1)
summary(mdl1)

mdl1_noC <- lm(SSP1$Lchl[c(2:6)] ~ Hist1$Lchl[c(2:6)])
ypred1_noC <- predict.lm(mdl1_noC)
summary(mdl1_noC)

mdl2 <- lm(SSP2$Lchl ~ Hist2$Lchl[1:6])
ypred2 <- predict.lm(mdl2)
summary(mdl2)

mdl2_noC <- lm(SSP2$Lchl[2:6] ~ Hist2$Lchl[2:6])
ypred2_noC <- predict.lm(mdl2_noC)
summary(mdl2_noC)

mdl3 <- lm(SSP3$Lchl ~ Hist3$Lchl[1:6])
ypred3 <- predict.lm(mdl3)
summary(mdl3)

mdl3_noC <- lm(SSP3$Lchl[2:6] ~ Hist3$Lchl[2:6])
ypred3_noC <- predict.lm(mdl3_noC)
summary(mdl3_noC)

Cff <- as.data.frame(mdl1$coefficients)
Cff[,2] <- mdl2$coefficients
Cff[,3] <- mdl3$coefficients
Cff[,4] <- mdl0$coefficients
Cff[,5] <- mdl1_noC$coefficients
Cff[,6] <- mdl2_noC$coefficients
Cff[,7] <- mdl3_noC$coefficients
Cff[,8] <- mdl0_noC$coefficients
names(Cff) <- c("LC","HCSS","HCPS","Global","LC_noC","HCSS_noC","HCPS_noC","Global_noC")

EC0 <- as.data.frame(Hist0$Lchl[1:6])
names(EC0) <- "Hist"
EC0$SSP <- SSP0$Lchl
EC0$model <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

EC0_noC <- EC0[c(2:6),]

EC1 <- as.data.frame(Hist1$Lchl[1:6])
names(EC1) <- "Hist"
EC1$SSP <- SSP1$Lchl
EC1$model <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

EC1_noC <- EC1[c(2:6),]

EC2 <- as.data.frame(Hist2$Lchl[1:6])
names(EC2) <- "Hist"
EC2$SSP <- SSP2$Lchl
EC2$model <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

EC2_noC <- EC2[c(2:6),]

EC3 <- as.data.frame(Hist3$Lchl[1:6])
names(EC3) <- "Hist"
EC3$SSP <- SSP3$Lchl
EC3$model <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

EC3_noC <- EC3[c(2:6),]

## EC plots
# Add 1:1 line!!!

c0 <- ggplot(EC0, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl Globally") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred0)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist0$Lchl[7], color = "grey") 

c1 <- ggplot(EC1, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in LC") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred1)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist1$Lchl[7], color = "grey") #+ 
  
c2 <- ggplot(EC2, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCSS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred2)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist2$Lchl[7], color = "grey")

c3 <- ggplot(EC3, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCPS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#228833","#999933","#33BBEE","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred3)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist3$Lchl[7], color = "grey")

c01 <- ggplot(EC0_noC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl Globally") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#33BBEE","#999933","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred0_noC)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist0$Lchl[7], color = "grey")

c11 <- ggplot(EC1_noC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in LC") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#33BBEE","#999933","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred1_noC)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist1$Lchl[7], color = "grey")

c21 <- ggplot(EC2_noC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCSS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#33BBEE","#999933","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred2_noC)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist2$Lchl[7], color = "grey")

c31 <- ggplot(EC3_noC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCPS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("#33BBEE","#999933","#004488","#EE6677","#AA3377")) +
  geom_line(aes(x = Hist, y = ypred3_noC)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist3$Lchl[7], color = "grey")

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_global_obsglm.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h0,s0,c0,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome1_LC45_obsglm.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h1,s1,c1,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome2_HCSS_obsglm.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h2,s2,c2,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome3_HCPS_obsglm.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h3,s3,c3,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

#without CAN outlier
pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_global_obsglm_noCAN.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h0,s0,c01,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome1_LC45_obsglm_noCAN.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h1,s1,c11,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome2_HCSS_obsglm_noCAN.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h2,s2,c21,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome3_HCPS_obsglm_noCAN.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( h3,s3,c31,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()




# ==========================================================================
# SST
# ==========================================================================

# load data
Hist <- read.csv("Hist_coeffs_mod_sst_trop_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
names(Hist) <- c("int","slop")

SSP <- read.csv("SSP585_coeffs_mod_sst_trop_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
names(SSP) <- c("int","slop")

## predict mesozoo over sst values
sst <- -2:30
hzoo <- matrix(NA, nrow = 5, ncol = 33)
for(i in 1:5)
{
  hzoo[i,] <- Hist$int[i] + Hist$slop[i] * sst
}

szoo <- matrix(NA, nrow = 5, ncol = 33)
for(i in 1:5)
{
  szoo[i,] <- SSP$int[i] + SSP$slop[i] * sst
}

hzoo <- as.data.frame(t(hzoo))
names(hzoo) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")
szoo <- as.data.frame(t(szoo))
names(szoo) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

hzoo$sst <- sst
szoo$sst <- sst


Hdf <- hzoo %>%
  select(sst,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -sst)

Sdf <- szoo %>%
  select(sst,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -sst)

## figures
b1 <- ggplot(Hdf, aes(x = sst, y = value)) + 
  geom_line(aes(color = model)) + 
  ggtitle("Hist") + 
  ylab("zmeso gC m-2") + xlab("sst C") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

b2 <- ggplot(Sdf, aes(x = sst, y = value)) + 
  geom_line(aes(color = model)) + 
  ggtitle("SSP 585") + 
  ylab("zmeso gC m-2") + xlab("sst C") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

pdf(file = 'Hist_SSP585_lm_zm_sst_trop_all_clim_200.pdf', width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( b1,b2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()



### linear regression of hist vs. future
mdl <- lm(SSP$slop ~ Hist$slop)
ypred <- predict.lm(mdl)
summary(mdl)
#Multiple R-squared:  0.8489,	Adjusted R-squared:  0.7985 
#F-statistic: 16.85 on 1 and 3 DF,  p-value: 0.02617

Cff <- as.data.frame(mdl$coefficients)
names(Cff) <- "all"

EC <- as.data.frame(Hist$slop)
names(EC) <- "Hist"
EC$SSP <- SSP$slop
EC$model <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM")

##

c1 <- ggplot(EC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to sst") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred)) + 
  annotate( geom = 'text', y = 20, x = -25, hjust = 0, label="r = 0.80", size=5) +
  annotate( geom = 'text', y = 18, x = -25, hjust = 0, label="p = 0.03", size=5)


pdf(file = 'Hist_SSP585_lm_zm_sst_EC_trop_all_clim_200.pdf', width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( b1,b2,c1,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


