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

# load data
Hist1 <- read.csv("Hist_coeffs_mod_chl_biome1_LC_obsglm.csv",sep=",",header = T,stringsAsFactors = F)
Hist2 <- read.csv("Hist_coeffs_mod_chl_biome2_HCSS_obsglm.csv",sep=",",header = T,stringsAsFactors = F)
Hist3 <- read.csv("Hist_coeffs_mod_chl_biome3_HCPS_obsglm.csv",sep=",",header = T,stringsAsFactors = F)
Herr <- read.csv("Hist_stderr_coeff_obsglm_chl_all_biomes.csv",sep=",",header = T,stringsAsFactors = F)
Oerr1 <- read.csv("Hist_stderr_predict_obsglm_chl_biome1.csv",sep=",",header = T,stringsAsFactors = F)
Oerr2 <- read.csv("Hist_stderr_predict_obsglm_chl_biome2.csv",sep=",",header = T,stringsAsFactors = F)
Oerr3 <- read.csv("Hist_stderr_predict_obsglm_chl_biome3.csv",sep=",",header = T,stringsAsFactors = F)

SSP1 <- read.csv("SSP585_coeffs_mod_chl_biome1_LC.csv",sep=",",header = T,stringsAsFactors = F)
SSP2 <- read.csv("SSP585_coeffs_mod_chl_biome2_HCSS.csv",sep=",",header = T,stringsAsFactors = F)
SSP3 <- read.csv("SSP585_coeffs_mod_chl_biome3_HCPS.csv",sep=",",header = T,stringsAsFactors = F)

## predict mesozoo over chl values
chl <- logspace(-5,-2,20)
#hzoo <- matrix(NA, nrow = 6, ncol = 20) #6 with obs
hzoo1 <- matrix(NA, nrow = 6, ncol = 20)
hzoo2 <- matrix(NA, nrow = 6, ncol = 20)
hzoo3 <- matrix(NA, nrow = 6, ncol = 20)
for(i in 1:6)
{
  hzoo1[i,] <- Hist1$Intercept[i] + Hist1$Lchl[i] * log(chl)
  hzoo2[i,] <- Hist2$Intercept[i] + Hist2$Lchl[i] * log(chl)
  hzoo3[i,] <- Hist3$Intercept[i] + Hist3$Lchl[i] * log(chl)
}

szoo1 <- matrix(NA, nrow = 5, ncol = 20)
szoo2 <- matrix(NA, nrow = 5, ncol = 20)
szoo3 <- matrix(NA, nrow = 5, ncol = 20)
for(i in 1:5)
{
  szoo1[i,] <- SSP1$Intercept[i] + SSP1$Lchl[i] * log(chl)
  szoo2[i,] <- SSP2$Intercept[i] + SSP2$Lchl[i] * log(chl)
  szoo3[i,] <- SSP3$Intercept[i] + SSP3$Lchl[i] * log(chl)
}

hzoo1 <- as.data.frame(t(hzoo1))
hzoo2 <- as.data.frame(t(hzoo2))
hzoo3 <- as.data.frame(t(hzoo3))
#names(hzoo) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo1) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo2) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
names(hzoo3) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")

szoo1 <- as.data.frame(t(szoo1))
szoo2 <- as.data.frame(t(szoo2))
szoo3 <- as.data.frame(t(szoo3))
names(szoo1) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
names(szoo2) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
names(szoo3) <- c("CAN","CNRM","GFDL","IPSL","UKESM")

hzoo1$chl <- chl
hzoo2$chl <- chl
hzoo3$chl <- chl
szoo1$chl <- chl
szoo2$chl <- chl
szoo3$chl <- chl

# Add std error high and low
hzoo1$high <- Hist1$Intercept[6]+Herr$Intercept[1] + (Hist1$Lchl[6]+Herr$Lchl[1]) * log(chl)
hzoo1$low <- Hist1$Intercept[6]-Herr$Intercept[1] + (Hist1$Lchl[6]-Herr$Lchl[1]) * log(chl)
hzoo2$high <- Hist2$Intercept[6]+Herr$Intercept[2] + (Hist2$Lchl[6]+Herr$Lchl[2]) * log(chl)
hzoo2$low <- Hist2$Intercept[6]-Herr$Intercept[2] + (Hist2$Lchl[6]-Herr$Lchl[2]) * log(chl)
hzoo3$high <- Hist3$Intercept[6]+Herr$Intercept[3] + (Hist3$Lchl[6]+Herr$Lchl[3]) * log(chl)
hzoo3$low <- Hist3$Intercept[6]-Herr$Intercept[3] + (Hist3$Lchl[6]-Herr$Lchl[2]) * log(chl)

# hzoo1$upr <- Oerr1$upr
# hzoo1$lwr <- Oerr1$lwr
# hzoo2$upr <- Oerr2$upr
# hzoo2$lwr <- Oerr2$lwr
# hzoo3$upr <- Oerr3$upr
# hzoo3$lwr <- Oerr3$lwr

# Hdf <- hzoo %>%
#   select(chl,CAN,CNRM,GFDL,IPSL,UKESM,obs) %>%
#   gather(key = "model", value = "value", -chl)
Hdf1 <- hzoo1 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM,obs,high,low) %>%
  gather(key = "model", value = "value", -chl)
Hdf2 <- hzoo2 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM,obs,high,low) %>%
  gather(key = "model", value = "value", -chl)
Hdf3 <- hzoo3 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM,obs,high,low) %>%
  gather(key = "model", value = "value", -chl)

Sdf1 <- szoo1 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)
Sdf2 <- szoo2 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)
Sdf3 <- szoo3 %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)

## figures
ylmts <- c(-12,8)
h1 <- ggplot(Hdf1, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + ylim(ylmts) +
  scale_x_log10() + ggtitle("LC Hist") + 
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","grey","red","grey","black","maroon"))
  #scale_color_manual(values = c("green4","blue","purple","red","black","maroon"))

s1 <- ggplot(Sdf1, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + ylim(ylmts) +
  scale_x_log10() + ggtitle("LC SSP 585") + 
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_biome1_LC_obsglm.pdf'), width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( h1,s1,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

ylmts2 <- c(-10,5)
h2 <- ggplot(Hdf2, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCSS Hist") + ylim(ylmts2) +
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","grey","red","grey","black","maroon"))
  #scale_color_manual(values = c("green4","blue","purple","red","maroon"))

s2 <- ggplot(Sdf2, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCSS SSP 585") + ylim(ylmts2) +
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_biome2_HCSS_obsglm.pdf'), width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( h2,s2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

h3 <- ggplot(Hdf3, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCPS Hist") + ylim(ylmts2) +
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","grey","red","grey","black","maroon"))
  #scale_color_manual(values = c("green4","blue","purple","red","maroon"))

s3 <- ggplot(Sdf3, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("HCPS SSP 585") + ylim(ylmts2) +
  ylab("zmeso gC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_biome3_HCPS_obsglm.pdf'), width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( h3,s3,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()



### linear regression of hist vs. future -----------------------------
mdl1 <- lm(SSP1$Lchl ~ Hist1$Lchl[1:5])
ypred1 <- predict.lm(mdl1)
summary(mdl1)

mdl1_noI <- lm(SSP1$Lchl[c(1:3,5)] ~ Hist1$Lchl[c(1:3,5)])
ypred1_noI <- predict.lm(mdl1_noI)
summary(mdl1_noI)

mdl2 <- lm(SSP2$Lchl ~ Hist2$Lchl[1:5])
ypred2 <- predict.lm(mdl2)
summary(mdl2)

mdl2_noC <- lm(SSP2$Lchl[2:5] ~ Hist2$Lchl[2:5])
ypred2_noC <- predict.lm(mdl2_noC)
summary(mdl2_noC)

mdl3 <- lm(SSP3$Lchl ~ Hist3$Lchl[1:5])
ypred3 <- predict.lm(mdl3)
summary(mdl3)

mdl3_noU <- lm(SSP3$Lchl[1:4] ~ Hist3$Lchl[1:4])
ypred3_noU <- predict.lm(mdl3_noU)
summary(mdl3_noU)

Cff <- as.data.frame(mdl1$coefficients)
Cff[,2] <- mdl2$coefficients
Cff[,3] <- mdl3$coefficients
Cff[,4] <- mdl1_noI$coefficients
Cff[,5] <- mdl2_noC$coefficients
Cff[,6] <- mdl3_noU$coefficients
names(Cff) <- c("LC","HCSS","HCPS","LC_noI","HCSS_noC","HCPS_noU")

EC1 <- as.data.frame(Hist1$Lchl[1:5])
names(EC1) <- "Hist"
EC1$SSP <- SSP1$Lchl
EC1$model <- c("CAN","CNRM","GFDL","IPSL","UKESM")

EC1_noI <- EC1[c(1:3,5),]

EC2 <- as.data.frame(Hist2$Lchl[1:5])
names(EC2) <- "Hist"
EC2$SSP <- SSP2$Lchl
EC2$model <- c("CAN","CNRM","GFDL","IPSL","UKESM")

EC2_noC <- EC2[2:5,]

EC3 <- as.data.frame(Hist3$Lchl[1:5])
names(EC3) <- "Hist"
EC3$SSP <- SSP3$Lchl
EC3$model <- c("CAN","CNRM","GFDL","IPSL","UKESM")

EC3_noU <- EC3[1:4,]

## EC plots
# Add 1:1 line!!!

c1 <- ggplot(EC1, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in LC") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred1)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist1$Lchl[6], color = "grey") #+ 
  # annotate( geom = 'text', y = 10, x = 1, hjust = 0, label="r = 0.99", size=5) +
  # annotate( geom = 'text', y = 9, x = 1, hjust = 0, label="p = 0.0002", size=5)

c2 <- ggplot(EC2, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCSS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred2)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist2$Lchl[6], color = "grey")

c3 <- ggplot(EC3, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCPS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred3)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist3$Lchl[6], color = "grey")

c11 <- ggplot(EC1_noI, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in LC") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","maroon")) +
  geom_line(aes(x = Hist, y = ypred1_noI)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist1$Lchl[6], color = "grey")

c21 <- ggplot(EC2_noC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCSS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred2_noC)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist2$Lchl[6], color = "grey")

c31 <- ggplot(EC3_noU, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl in HCPS") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red")) +
  geom_line(aes(x = Hist, y = ypred3_noU)) +
  geom_abline(intercept = 0, linetype="dashed", color = "black") +
  geom_vline(xintercept = Hist3$Lchl[6], color = "grey")

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome1_LC_obsglm.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
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

#without outliers
pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome1_LC_obsglm_noIPSL.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
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

pdf(file = paste0(figp,'Hist_SSP585_lm_zm_chl_EC_biome3_HCPS_obsglm_noUK.pdf'), width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
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
names(hzoo) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
szoo <- as.data.frame(t(szoo))
names(szoo) <- c("CAN","CNRM","GFDL","IPSL","UKESM")

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
EC$model <- c("CAN","CNRM","GFDL","IPSL","UKESM")

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


