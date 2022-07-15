# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes
# GLMM100

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
Zm <- read.csv(paste0(datap,"skill_hist_model_obsglm100_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv(paste0(datap,"hist_chl_model_obsglm100_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)
Sst <- read.csv(paste0(datap,"hist_sst_model_obsglm100_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)

SZm <- read.csv(paste0(datap,"skill_hist_Stromberg_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
SChl <- read.csv(paste0(datap,"hist_chl_seawifs_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)


CM <- as.data.frame(Zm[,c(1:3,10)])
CM[,5] <- as.data.frame(Chl[,3])
CM[,6] <- as.data.frame(Sst[,3])
names(CM) <- c("Lat","Lon","zmeso","biome","chl","sst")
CM <- na.omit(CM)

MM <- as.data.frame(Zm[,c(1:2,4,11)])
MM[,5] <- as.data.frame(Chl[,4])
MM[,6] <- as.data.frame(Sst[,4])
names(MM) <- c("Lat","Lon","zmeso","biome","chl","sst")
MM <- na.omit(MM)

NM <- as.data.frame(Zm[,c(1:2,5,12)])
NM[,5] <- as.data.frame(Chl[,5])
NM[,6] <- as.data.frame(Sst[,5])
names(NM) <- c("Lat","Lon","zmeso","biome","chl","sst")
NM <- na.omit(NM)

GM <- as.data.frame(Zm[,c(1:2,6,13)])
GM[,5] <- as.data.frame(Chl[,6])
GM[,6] <- as.data.frame(Sst[,6])
names(GM) <- c("Lat","Lon","zmeso","biome","chl","sst")
GM <- na.omit(GM)

IM <- as.data.frame(Zm[,c(1:2,7,14)])
IM[,5] <- as.data.frame(Chl[,7])
IM[,6] <- as.data.frame(Sst[,7])
names(IM) <- c("Lat","Lon","zmeso","biome","chl","sst")
IM <- na.omit(IM)

UM <- as.data.frame(Zm[,c(1:2,8,15)])
UM[,5] <- as.data.frame(Chl[,8])
UM[,6] <- as.data.frame(Sst[,8])
names(UM) <- c("Lat","Lon","zmeso","biome","chl","sst")
UM <- na.omit(UM)

OG <- as.data.frame(Zm[,c(1:2,9,16)])
OG[,5] <- as.data.frame(Chl[,9])
OG[,6] <- as.data.frame(Sst[,9])
names(OG) <- c("Lat","Lon","zmeso","biome","chl","sst")
OG <- na.omit(OG)

SM <- as.data.frame(SZm)
SM[,5] <- as.data.frame(SChl[,3])
names(SM) <- c("Lat","Lon","zmeso","biome","chl")
SM <- na.omit(SM)

CM$Lzmeso <- log10(CM$zmeso+1e-6)
MM$Lzmeso <- log10(MM$zmeso+1e-6)
NM$Lzmeso <- log10(NM$zmeso+1e-6)
GM$Lzmeso <- log10(GM$zmeso+1e-6)
IM$Lzmeso <- log10(IM$zmeso+1e-6)
UM$Lzmeso <- log10(UM$zmeso+1e-6)
OG$Lzmeso <- log10(OG$zmeso+1e-6)
SM$Lzmeso <- log10(SM$zmeso+1e-6)

CM$Lchl <- log10(CM$chl+1e-6)
MM$Lchl <- log10(MM$chl+1e-6)
NM$Lchl <- log10(NM$chl+1e-6)
GM$Lchl <- log10(GM$chl+1e-6)
IM$Lchl <- log10(IM$chl+1e-6)
UM$Lchl <- log10(UM$chl+1e-6)
OG$Lchl <- log10(OG$chl+1e-6)
SM$Lchl <- log10(SM$chl+1e-6)

## Restrict latitudes to those observed
maxlat <- max(OG$Lat)
minlat <- min(OG$Lat)
CM <- subset(CM, Lat >= minlat & Lat <= maxlat)
MM <- subset(MM, Lat >= minlat & Lat <= maxlat)
NM <- subset(NM, Lat >= minlat & Lat <= maxlat)
GM <- subset(GM, Lat >= minlat & Lat <= maxlat)
IM <- subset(IM, Lat >= minlat & Lat <= maxlat)
UM <- subset(UM, Lat >= minlat & Lat <= maxlat)

## Remove LC >45 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
cid1 <- which(CM$biome==1)
cid2 <- which(CM$Lat > 45 | CM$Lat < -45)
cid <- intersect(cid1,cid2)
CM <- CM[-cid,]

mid1 <- which(MM$biome==1)
mid2 <- which(MM$Lat > 45 | MM$Lat < -45)
mid <- intersect(mid1,mid2)
MM <- MM[-mid,]

nid1 <- which(NM$biome==1)
nid2 <- which(NM$Lat > 45 | NM$Lat < -45)
nid <- intersect(nid1,nid2)
NM <- NM[-nid,]

gid1 <- which(GM$biome==1)
gid2 <- which(GM$Lat > 45 | GM$Lat < -45)
gid <- intersect(gid1,gid2)
GM <- GM[-gid,]

iid1 <- which(IM$biome==1)
iid2 <- which(IM$Lat > 45 | IM$Lat < -45)
iid <- intersect(iid1,iid2)
IM <- IM[-iid,]

uid1 <- which(UM$biome==1)
uid2 <- which(UM$Lat > 45 | UM$Lat < -45)
uid <- intersect(uid1,uid2)
UM <- UM[-uid,]

oid1 <- which(OG$biome==1)
oid2 <- which(OG$Lat > 45 | OG$Lat < -45)
oid <- intersect(oid1,oid2)
OG <- OG[-oid,]

sid1 <- which(SM$biome==1)
sid2 <- which(SM$Lat > 45 | SM$Lat < -45)
sid <- intersect(sid1,sid2)
SM <- SM[-sid,]


### Plots -------------------------------------------------------
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CAN zmeso") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c2 <- ggplot(subset(CM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c3 <- ggplot(subset(CM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c0 <- ggplot(CM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

m1 <- ggplot(subset(MM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CMCC zmeso") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10(limits=c(0.02,0.15)) + theme(legend.position='none')
m2 <- ggplot(subset(MM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m3 <- ggplot(subset(MM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m0 <- ggplot(MM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 500) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

n1 <- ggplot(subset(NM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CNRM zmeso") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n2 <- ggplot(subset(NM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n3 <- ggplot(subset(NM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n0 <- ggplot(NM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
  
g1 <- ggplot(subset(GM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("GFDL zmeso") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g2 <- ggplot(subset(GM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g3 <- ggplot(subset(GM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g0 <- ggplot(GM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

i1 <- ggplot(subset(IM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("IPSL zmeso") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i2 <- ggplot(subset(IM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i3 <- ggplot(subset(IM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i0 <- ggplot(IM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

u1 <- ggplot(subset(UM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + #, se=FALSE, col="blue", size = 0.25) +  
  ylab("UK zmeso") + xlab("chl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u2 <- ggplot(subset(UM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u3 <- ggplot(subset(UM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u0 <- ggplot(UM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

o1 <- ggplot(subset(OG,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("obsGLMM") + xlab("Mchl") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o2 <- ggplot(subset(OG,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o3 <- ggplot(subset(OG,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o0 <- ggplot(OG, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

s1 <- ggplot(subset(SM,biome==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + #, se=FALSE, col="blue", size = 0.25) +  
  ylab("obsSM") + xlab("Schl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s2 <- ggplot(subset(SM,biome==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s3 <- ggplot(subset(SM,biome==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s0 <- ggplot(SM, aes(y=zmeso, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')


# COPEPOD OBS
Zm <- read.csv(paste0(datap,"skill_hist_COPEPOD_all_clim_200_chls.csv"),sep=",",header = T,stringsAsFactors = F)

CSM <- as.data.frame(Zm[,c("Lat","Lon","zmeso200","SEAWIFSchl","SEAWIFSbiomes")])
CSM <- na.omit(CSM)
names(CSM) <- c("Lat","Lon","zmeso200","chl","biome")

CMM <- as.data.frame(Zm[,c("Lat","Lon","zmeso200","MODISchl","MODISbiomes")])
CMM <- na.omit(CMM)
names(CMM) <- c("Lat","Lon","zmeso200","chl","biome")

# from mgC to gC
CSM$zmeso200 <- CSM$zmeso200*1e-3
CMM$zmeso200 <- CMM$zmeso200*1e-3

CSM$Lzmeso <- log10(CSM$zmeso200+1e-6)
CMM$Lzmeso <- log10(CMM$zmeso200+1e-6)

CSM$Lchl <- log10(CSM$chl+1e-6)
CMM$Lchl <- log10(CMM$chl+1e-6)

# Remove LC >45 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
sid1 <- which(CSM$biome==1)
sid2 <- which(CSM$Lat > 45 | CSM$Lat < -45)
sid <- intersect(sid1,sid2)
CSM <- CSM[-sid,]

mid1 <- which(CMM$biome==1)
mid2 <- which(CMM$Lat > 45 | CMM$Lat < -45)
mid <- intersect(mid1,mid2)
CMM <- CMM[-mid,]

cm1 <- ggplot(subset(CMM,biome==1), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("COPEPOD") + xlab("Mchl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10(limits=c(0.02,0.15)) + theme(legend.position='none')
cm2 <- ggplot(subset(CMM,biome==2), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
cm3 <- ggplot(subset(CMM,biome==3), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
cm0 <- ggplot(CMM, aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("Mchl") + #ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

cs1 <- ggplot(subset(CSM,biome==1), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + #, se=FALSE, col="blue", size = 0.25) +  
  ylab("COPEPOD") + xlab("Schl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
cs2 <- ggplot(subset(CSM,biome==2), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
cs3 <- ggplot(subset(CSM,biome==3), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
cs0 <- ggplot(CSM, aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("Schl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')


## With Global
png(paste0(figp,'scatter_hist_zmeso_chl_biomes_global_all_log_LC45.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m0,
           NULL,NULL,NULL,NULL,
           n1,n2,n3,n0,
           NULL,NULL,NULL,NULL,
           g1,g2,g3,g0,
           NULL,NULL,NULL,NULL,
           i1,i2,i3,i0,
           NULL,NULL,NULL,NULL,
           u1,u2,u3,u0,
           nrow = 9, ncol = 4,
           rel_widths = c( 1, 1, 1, 1 ), 
           rel_heights = c(1,-0.15,1,-0.15,1,-0.15,1,-0.15,1) ,
           align = 'h' )
dev.off()

png(paste0(figp,'scatter_zmeso_chl_biomes_global_all_log_glm100_strom_copepod_LC45.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( o1,o2,o3,o0,
           NULL,NULL,NULL,NULL,
           s1,s2,s3,s0,
           NULL,NULL,NULL,NULL,
           cm1,cm2,cm3,cm0,
           NULL,NULL,NULL,NULL,
           cs1,cs2,cs3,cs0,
           nrow = 7, ncol = 4,
           rel_widths = c( 1, 1, 1, 1 ), 
           rel_heights = c(1,-0.05,1,-0.05,1,-0.05,1) ,
           align = 'h' )
dev.off()


