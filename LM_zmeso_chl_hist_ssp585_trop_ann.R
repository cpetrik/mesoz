# Linear regressions of zmeso biomass with surf chl 
# Tropical regions, annual climatology
# Historic and SSP585

rm(list=ls())

library(plyr)
library(cowplot) #plot_grid
library(R.matlab)
library(rdetools)
library("tidyverse")

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
Hist <- read.csv("Hist_coeffs_mod_chl_trop_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
names(Hist)[1] <- "model"

SSP <- read.csv("SSP585_coeffs_mod_chl_trop_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
names(SSP)[1] <- "model"

## predict mesozoo over chl values
chl <- logspace(-5,-2,20)
hzoo <- matrix(NA, nrow = 6, ncol = 20)
for(i in 1:6)
{
  hzoo[i,] <- Hist$int[i] + Hist$slop[i] * log10(chl)
}

szoo <- matrix(NA, nrow = 5, ncol = 20)
for(i in 1:5)
{
  szoo[i,] <- SSP$int[i] + SSP$slop[i] * log10(chl)
}

hzoo <- as.data.frame(t(hzoo))
names(hzoo) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
szoo <- as.data.frame(t(szoo))
names(szoo) <- c("CAN","CNRM","GFDL","IPSL","UKESM")

hzoo$chl <- chl
szoo$chl <- chl

Hdf <- hzoo %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM,obs) %>%
  gather(key = "model", value = "value", -chl)

Sdf <- szoo %>%
  select(chl,CAN,CNRM,GFDL,IPSL,UKESM) %>%
  gather(key = "model", value = "value", -chl)

## figures
b1 <- ggplot(Hdf, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("Hist") + 
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","red","black","maroon"))

b2 <- ggplot(Sdf, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("SSP 585") + 
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

a1 <- ggplot(Hdf, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("Hist") + 
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  ylim(0,5) +
  scale_color_manual(values = c("green4","blue","purple","red","black","maroon"))

a2 <- ggplot(Sdf, aes(x = chl, y = value)) + 
  geom_line(aes(color = model)) + 
  scale_x_log10() + ggtitle("SSP 585") + 
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  ylim(0,5) + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

pdf(file = 'Hist_SSP585_lm_zm_chl_trop_all_clim_200.pdf', width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( b1,b2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = 'Hist_SSP585_lm_zm_chl_trop_all_clim_200_noCAN.pdf', width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( a1,a2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### linear regression of hist vs. future
mdl <- lm(SSP$slop ~ Hist$slop[1:5])
ypred <- predict.lm(mdl)
summary(mdl)
#Multiple R-squared:  0.9926,	Adjusted R-squared:  0.9901 
#F-statistic: 399.8 on 1 and 3 DF,  p-value: 0.0002734

mdl2 <- lm(SSP$slop[2:5] ~ Hist$slop[2:5])
ypred2 <- predict(mdl2)
summary(mdl2)
#Multiple R-squared:  0.3633,	Adjusted R-squared:  0.04493 
#F-statistic: 1.141 on 1 and 2 DF,  p-value: 0.3973

Cff <- as.data.frame(mdl$coefficients)
Cff[,2] <- mdl2$coefficients
names(Cff) <- c("all","noCAN")

EC <- as.data.frame(Hist$slop[1:5])
names(EC) <- "Hist"
EC$SSP <- SSP$slop
EC$model <- c("CAN","CNRM","GFDL","IPSL","UKESM")

##

c1 <- ggplot(EC, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred)) + 
  annotate( geom = 'text', y = 10, x = 1, hjust = 0, label="r = 0.99", size=5) +
  annotate( geom = 'text', y = 9, x = 1, hjust = 0, label="p = 0.0002", size=5)


EC2 <- EC[2:5,]
c2 <- ggplot(EC2, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred2))+ 
  annotate( geom = 'text', y = 1.5, x = 0.6, hjust = 0, label="r = 0.04", size=5) +
  annotate( geom = 'text', y = 1.4, x = 0.6, hjust = 0, label="p = 0.40", size=5)


pdf(file = 'Hist_SSP585_lm_zm_chl_EC_trop_all_clim_200.pdf', width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( b1,b2,c1,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

pdf(file = 'Hist_SSP585_lm_zm_chl_EC_trop_all_clim_200_noCAN.pdf', width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( a1,a2,c2,
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
  ylab("zmeso mgC m-2") + xlab("sst C") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))

b2 <- ggplot(Sdf, aes(x = sst, y = value)) + 
  geom_line(aes(color = model)) + 
  ggtitle("SSP 585") + 
  ylab("zmeso mgC m-2") + xlab("sst C") + 
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


