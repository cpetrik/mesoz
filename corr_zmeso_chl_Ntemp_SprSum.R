# Calculate correlation of modeled zmeso biomass and chl

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
CM <- read.csv("can_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Ann <- read.csv("obs_mod_chl_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## Latitude < 60 N & >= 30 N: JJA + MAM
aCM <- subset(CM, Lat > 30 & Lat < 60)
aNM <- subset(NM, Lat > 30 & Lat < 60)
aGM <- subset(GM, Lat > 30 & Lat < 60)
aIM <- subset(IM, Lat > 30 & Lat < 60)
aUM <- subset(UM, Lat > 30 & Lat < 60)
Amz <- subset(Ann, Lat > 30 & Lat < 60)

## Take mean over MAM and JJA
aCM$chl <- (aCM$chlMAM+aCM$chlJJA)/2
aNM$chl <- (aNM$chlMAM+aNM$chlJJA)/2
aGM$chl <- (aGM$chlMAM+aGM$chlJJA)/2
aIM$chl <- (aIM$chlMAM+aIM$chlJJA)/2
aUM$chl <- (aUM$chlMAM+aUM$chlJJA)/2

aCM$sst <- (aCM$sstMAM+aCM$sstJJA)/2
aNM$sst <- (aNM$sstMAM+aNM$sstJJA)/2
aGM$sst <- (aGM$sstMAM+aGM$sstJJA)/2
aIM$sst <- (aIM$sstMAM+aIM$sstJJA)/2
aUM$sst <- (aUM$sstMAM+aUM$sstJJA)/2

aCM$zoo <- (aCM$zooMAM+aCM$zooJJA)/2
aNM$zoo <- (aNM$zooMAM+aNM$zooJJA)/2
aGM$zoo <- (aGM$zooMAM+aGM$zooJJA)/2
aIM$zoo <- (aIM$zooMAM+aIM$zooJJA)/2
aUM$zoo <- (aUM$zooMAM+aUM$zooJJA)/2

## Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chl * 1e3
aGM$chl <- aGM$chl * 1e3
aUM$chl <- aUM$chl * 1e3

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,c("zoo","chl","sst")])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
iid <- !is.infinite(rowSums(Cmydata))
Cmydata <- Cmydata[iid,]
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,c("zoo","chl","sst")])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,c("zoo","chl","sst")])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,c("zoo","chl","sst")])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,c("zoo","chl","sst")])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)

Amydata <- as.matrix(Amz[,c("obs","chl")])
Amydata <- na.omit(Amydata)
Acorr <- rcorr(Amydata) #,type="pearson","spearman"
Acorr_all_p <- data.frame(Acorr$P)
Acorr_all_r <- data.frame(Acorr$r)

### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl
#get rid of infinite when trans
aCM <- aCM[iid,]

a1 <- ggplot(aCM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
   annotate( geom = 'text', y = 1e3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
   annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

# untransform observations
Amz$Robs <- 10^Amz$obs
Amz$Rchl <- 10^Amz$chl * 1e-3 #from mg to g

a6 <- ggplot(Amz, aes(y=Robs, x=Rchl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("Obs zmeso mgC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_mesoz_chl_Ntemp_SprSum_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,a6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  annotate( geom = 'text', y = 600, x = 10, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 550, x = 10, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5)
# scale_y_log10() + scale_x_log10() + 
# annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
# annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

s2 <- ggplot(aNM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  annotate( geom = 'text', y = 6000, x = 10, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 10, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  annotate( geom = 'text', y = 4000, x = 10, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3500, x = 10, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  annotate( geom = 'text', y = 4000, x = 10, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3500, x = 10, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  annotate( geom = 'text', y = 6000, x = 10, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 10, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_mesoz_sst_Ntemp_SprSum_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### LM of Chl ---------------------------------------------------------------
##Kwiatkowshki paper does not force through origin

# CAN
#check for outliers
boxplot(log10(aCM$zoo), main="Zmeso", 
        sub=paste("Outlier rows: ", boxplot.stats(log10(aCM$zoo))$out))
#distribution
library(e1071)
plot(density(log10(aCM$zoo)), main="Density Plot: Zmeso", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(log10(aCM$zoo)), 2)))  # density plot 

#correlation
cor(log10(aCM$zoo),log10(aCM$chl))

Clm <- lm(log10(zoo)~log10(chl), data=aCM)
summary(Clm)
Clm0 <- lm(log10(zoo)~log10(chl)-1, data=aCM)
summary(Clm0)
#predict.lm(Clm0,se.fit = T)




Nlm <- lm(log10(zoo)~log10(chl), data=aNM)
summary(Nlm)
Nlm0 <- lm(log10(zoo)~log10(chl)-1, data=aNM)
summary(Nlm0)

Glm <- lm(log10(zoo)~log10(chl), data=aGM)
summary(Glm)
Glm0 <- lm(log10(zoo)~log10(chl)-1, data=aGM)
summary(Glm0)

Ilm <- lm(log10(zoo)~log10(chl), data=aIM)
summary(Ilm)
Ilm0 <- lm(log10(zoo)~log10(chl)-1, data=aIM)
summary(Ilm0)

Ulm <- lm(log10(zoo)~log10(chl), data=aUM)
summary(Ulm)
Ulm0 <- lm(log10(zoo)~log10(chl)-1, data=aUM)
summary(Ulm0)


Olm <- lm(log10(Robs)~log10(Rchl), data=Amz)
summary(Olm)

#SST
ClmS <- lm(zoo~sst, data=aCM)
summary(ClmS)

NlmS <- lm(zoo~sst, data=aNM)
summary(NlmS)

GlmS <- lm(zoo~sst, data=aGM)
summary(GlmS)

IlmS <- lm(zoo~sst, data=aIM)
summary(IlmS)

UlmS <- lm(zoo~sst, data=aUM)
summary(UlmS)


## Save coefficients
Cff <- as.data.frame(t(coefficients(Clm)))
Cff[2,] <- t(coefficients(Nlm))
Cff[3,] <- t(coefficients(Glm))
Cff[4,] <- t(coefficients(Ilm))
Cff[5,] <- t(coefficients(Ulm))
Cff[6,] <- t(coefficients(Olm))
row.names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff,"Hist_coeffs_mod_chl_Ntemp_SprSum_200.csv",sep=",",row.names=T)

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"Hist_coeffs_mod_sst_Ntemp_SprSum_200.csv",sep=",",row.names=T)


# =========================================================================
#  SSP 585 
# =========================================================================

# load data
CM <- read.csv("can_ssp585_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_ssp585_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_ssp585_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_ssp585_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_ssp585_mod_zoo_chl_sst_SprSum_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## Latitude < 60 N & >= 30 N: JJA + MAM
aCM <- subset(CM, Lat > 30 & Lat < 60)
aNM <- subset(NM, Lat > 30 & Lat < 60)
aGM <- subset(GM, Lat > 30 & Lat < 60)
aIM <- subset(IM, Lat > 30 & Lat < 60)
aUM <- subset(UM, Lat > 30 & Lat < 60)

## Take mean over MAM and JJA
aCM$chl <- (aCM$chlMAM+aCM$chlJJA)/2
aNM$chl <- (aNM$chlMAM+aNM$chlJJA)/2
aGM$chl <- (aGM$chlMAM+aGM$chlJJA)/2
aIM$chl <- (aIM$chlMAM+aIM$chlJJA)/2
aUM$chl <- (aUM$chlMAM+aUM$chlJJA)/2

aCM$sst <- (aCM$sstMAM+aCM$sstJJA)/2
aNM$sst <- (aNM$sstMAM+aNM$sstJJA)/2
aGM$sst <- (aGM$sstMAM+aGM$sstJJA)/2
aIM$sst <- (aIM$sstMAM+aIM$sstJJA)/2
aUM$sst <- (aUM$sstMAM+aUM$sstJJA)/2

aCM$zoo <- (aCM$zooMAM+aCM$zooJJA)/2
aNM$zoo <- (aNM$zooMAM+aNM$zooJJA)/2
aGM$zoo <- (aGM$zooMAM+aGM$zooJJA)/2
aIM$zoo <- (aIM$zooMAM+aIM$zooJJA)/2
aUM$zoo <- (aUM$zooMAM+aUM$zooJJA)/2

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chl * 1e3
aGM$chl <- aGM$chl * 1e3
aUM$chl <- aUM$chl * 1e3

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,c("zoo","chl","sst")])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,c("zoo","chl","sst")])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,c("zoo","chl","sst")])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,c("zoo","chl","sst")])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,c("zoo","chl","sst")])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)



### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl
a1 <- ggplot(aCM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl kg m-3") + 
  # annotate( geom = 'text', y = 600, x = 0, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 550, x = 0, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e2, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 700, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zoo, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 2e-6, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e3, x = 2e-6, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_chl_Ntemp_SprSum_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  annotate( geom = 'text', y = 400, x = 5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 350, x = 5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5)
# scale_y_log10() + scale_x_log10() + 
# annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
# annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

s2 <- ggplot(aNM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3000, x = 0, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  annotate( geom = 'text', y = 3500, x = 17, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3300, x = 17, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  annotate( geom = 'text', y = 1000, x = 1.5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 900, x = 1.5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zoo, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  annotate( geom = 'text', y = 3500, x = 18, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3000, x = 18, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_sst_Ntemp_SprSum_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin, but are anomalies

# CAN
Clm <- lm(log10(zoo)~log10(chl), data=aCM)
summary(Clm)

Nlm <- lm(log10(zoo)~log10(chl), data=aNM)
summary(Nlm)

Glm <- lm(log10(zoo)~log10(chl), data=aGM)
summary(Glm)

Ilm <- lm(log10(zoo)~log10(chl), data=aIM)
summary(Ilm)

Ulm <- lm(log10(zoo)~log10(chl), data=aUM)
summary(Ulm)

#SST
ClmS <- lm(zoo~sst, data=aCM)
summary(ClmS)

NlmS <- lm(zoo~sst, data=aNM)
summary(NlmS)

GlmS <- lm(zoo~sst, data=aGM)
summary(GlmS)

IlmS <- lm(zoo~sst, data=aIM)
summary(IlmS)

UlmS <- lm(zoo~sst, data=aUM)
summary(UlmS)


## Save coefficients
Cff <- as.data.frame(t(coefficients(Clm)))
Cff[2,] <- t(coefficients(Nlm))
Cff[3,] <- t(coefficients(Glm))
Cff[4,] <- t(coefficients(Ilm))
Cff[5,] <- t(coefficients(Ulm))
row.names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Cff,"SSP585_coeffs_mod_chl_Ntemp_SprSum_200.csv",sep=",",row.names=T)

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"SSP585_coeffs_mod_sst_Ntemp_SprSum_200.csv",sep=",",row.names=T)


