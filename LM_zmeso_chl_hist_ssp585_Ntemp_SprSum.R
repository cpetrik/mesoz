# Linear regressions of zmeso biomass with surf chl 
# Ntemp regions, Spr-Sum
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
Hist <- read.csv("Hist_coeffs_mod_chl_Ntemp_SprSum_200.csv",sep=",",header = T,stringsAsFactors = F)
names(Hist)[1] <- "model"

SSP <- read.csv("SSP585_coeffs_mod_chl_Ntemp_SprSum_200.csv",sep=",",header = T,stringsAsFactors = F)
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
a1 <- ggplot(hzoo, aes( x=chl)) + theme_bw(base_size=14) +  
  geom_line(aes(y = CAN), color = "green4") + 
  geom_line(aes(y = CNRM), color="blue") + 
  geom_line(aes(y = GFDL), color = "purple") + 
  geom_line(aes(y = IPSL), color="red") + 
  geom_line(aes(y = UKESM), color="maroon") +
  geom_line(aes(y = obs), color="black") +
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  scale_x_log10() + ggtitle("Hist") + 
  theme(legend.text = names(hzoo)[1:6])

a2 <- ggplot(szoo, aes( x=chl)) + theme_bw(base_size=14) +  
  geom_line(aes(y = CAN), color = "green4") + 
  geom_line(aes(y = CNRM), color="blue") + 
  geom_line(aes(y = GFDL), color = "purple") + 
  geom_line(aes(y = IPSL), color="red") + 
  geom_line(aes(y = UKESM), color="maroon") + 
  ylab("zmeso mgC m-2") + xlab("chl g m-3") + 
  scale_x_log10() + ggtitle("SSP 585")


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

pdf(file = 'Hist_SSP585_lm_zm_chl_Ntemp_SprSum_clim_200.pdf', width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( b1,b2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### linear regression of hist vs. future
xvals <- 4:23
mdl <- lm(SSP$slop ~ Hist$slop[1:5])
ypred <- predict.lm(mdl)
summary(mdl)
#Multiple R-squared:  0.9865,	Adjusted R-squared:  0.982 
#F-statistic: 219.1 on 1 and 3 DF,  p-value: 0.0006691

mdl2 <- lm(SSP$slop[2:5] ~ Hist$slop[2:5])
ypred2 <- predict(mdl2)
summary(mdl2)
#Multiple R-squared:  0.9754,	Adjusted R-squared:  0.9631 
#F-statistic: 79.24 on 1 and 2 DF,  p-value: 0.01239

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
  annotate( geom = 'text', y = 7, x = 1, hjust = 0, label="r = 0.98", size=5) +
  annotate( geom = 'text', y = 6, x = 1, hjust = 0, label="p = 0.0007", size=5)


EC2 <- EC[2:5,]
c2 <- ggplot(EC2, aes(x = Hist, y = SSP)) + 
  geom_point(aes(color = model)) + 
  ggtitle("Climate sens of mesoz to chl") + 
  ylab("SSP 585") + xlab("Hist") + 
  scale_color_manual(values = c("blue","purple","red","maroon")) +
  geom_line(aes(x = Hist, y = ypred2))+ 
  annotate( geom = 'text', y = 2.5, x = 1, hjust = 0, label="r = 0.96", size=5) +
  annotate( geom = 'text', y = 2.3, x = 1, hjust = 0, label="p = 0.01", size=5)


pdf(file = 'Hist_SSP585_lm_zm_chl_EC_Ntemp_SprSum_clim_200.pdf', width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( b1,b2,c1,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


# ==========================================================================
# SST
# ==========================================================================

# load data
Hist <- read.csv("Hist_coeffs_mod_sst_Ntemp_SprSum_200.csv",sep=",",header = T,stringsAsFactors = F)
names(Hist) <- c("int","slop")

SSP <- read.csv("SSP585_coeffs_mod_sst_Ntemp_SprSum_200.csv",sep=",",header = T,stringsAsFactors = F)
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

pdf(file = 'Hist_SSP585_lm_zm_sst_Ntemp_SprSum_200.pdf', width = unit( 10, 'cm' ), height = unit( 5, 'cm' ))
plot_grid( b1,b2,
           nrow = 1, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()



### linear regression of hist vs. future
mdl <- lm(SSP$slop ~ Hist$slop)
ypred <- predict.lm(mdl)
summary(mdl)
#Multiple R-squared:  0.9682,	Adjusted R-squared:  0.9575 
#F-statistic:  91.2 on 1 and 3 DF,  p-value: 0.002436

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
  annotate( geom = 'text', y = 10, x = -50, hjust = 0, label="r = 0.96", size=5) +
  annotate( geom = 'text', y = 5, x = -50, hjust = 0, label="p = 0.002", size=5)


pdf(file = 'Hist_SSP585_lm_zm_sst_EC_Ntemp_SprSum_200.pdf', width = unit( 15, 'cm' ), height = unit( 4.5, 'cm' ))
plot_grid( b1,b2,c1,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

