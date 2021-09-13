# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Zm <- read.csv("ssp585_mesoz_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv("ssp585_chl_model_all_clim.csv",sep=",",header = T,stringsAsFactors = F)
Sst <- read.csv("ssp585_sst_model_all_clim.csv",sep=",",header = T,stringsAsFactors = F)

CM <- as.data.frame(Zm[,c(1:3,8)])
CM[,5] <- as.data.frame(Chl[,3])
CM[,6] <- as.data.frame(Sst[,3])
names(CM) <- c("Lat","Lon","mesoz","biome","chl","sst")
CM <- na.omit(CM)

NM <- as.data.frame(Zm[,c(1:2,4,9)])
NM[,5] <- as.data.frame(Chl[,4])
NM[,6] <- as.data.frame(Sst[,4])
names(NM) <- c("Lat","Lon","mesoz","biome","chl","sst")
NM <- na.omit(NM)

GM <- as.data.frame(Zm[,c(1:2,5,10)])
GM[,5] <- as.data.frame(Chl[,5])
GM[,6] <- as.data.frame(Sst[,5])
names(GM) <- c("Lat","Lon","mesoz","biome","chl","sst")
GM <- na.omit(GM)

IM <- as.data.frame(Zm[,c(1:2,6,11)])
IM[,5] <- as.data.frame(Chl[,6])
IM[,6] <- as.data.frame(Sst[,6])
names(IM) <- c("Lat","Lon","mesoz","biome","chl","sst")
IM <- na.omit(IM)

UM <- as.data.frame(Zm[,c(1:2,7,12)])
UM[,5] <- as.data.frame(Chl[,7])
UM[,6] <- as.data.frame(Sst[,7])
names(UM) <- c("Lat","Lon","mesoz","biome","chl","sst")
UM <- na.omit(UM)

CM$Lzmeso <- log(CM$mesoz+1e-16)
NM$Lzmeso <- log(NM$mesoz+1e-16)
GM$Lzmeso <- log(GM$mesoz+1e-16)
IM$Lzmeso <- log(IM$mesoz+1e-16)
UM$Lzmeso <- log(UM$mesoz+1e-16)

CM$Lchl <- log(CM$chl+1e-16)
NM$Lchl <- log(NM$chl+1e-16)
GM$Lchl <- log(GM$chl+1e-16)
IM$Lchl <- log(IM$chl+1e-16)
UM$Lchl <- log(UM$chl+1e-16)

## Restrict latitudes to those observed
oZm <- read.csv("skill_model_obsglm_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
oChl <- read.csv("chl_model_obsglm_all_clim.csv",sep=",",header = T,stringsAsFactors = F)
oSst <- read.csv("sst_model_obsglm_all_clim.csv",sep=",",header = T,stringsAsFactors = F)

OG <- as.data.frame(oZm[,c(1:2,8,14)])
OG[,5] <- as.data.frame(oChl[,8])
OG[,6] <- as.data.frame(oSst[,8])
names(OG) <- c("Lat","Lon","mesoz","biome","chl","sst")
OG <- na.omit(OG)

maxlat <- max(OG$Lat)
minlat <- min(OG$Lat)
CM <- subset(CM, Lat >= minlat & Lat <= maxlat)
NM <- subset(NM, Lat >= minlat & Lat <= maxlat)
GM <- subset(GM, Lat >= minlat & Lat <= maxlat)
IM <- subset(IM, Lat >= minlat & Lat <= maxlat)
UM <- subset(UM, Lat >= minlat & Lat <= maxlat)

## Remove LC >50 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
cid1 <- which(CM$biome==1)
cid2 <- which(CM$Lat > 45 | CM$Lat < -45)
cid <- intersect(cid1,cid2)
CM <- CM[-cid,]

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


### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Ccorr <- ddply(CM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ncorr <- ddply(NM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Gcorr <- ddply(GM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Icorr <- ddply(IM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ucorr <- ddply(UM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))


### Plots
## Zoo & Chl biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
c2 <- ggplot(subset(CM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
c3 <- ggplot(subset(CM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 

n1 <- ggplot(subset(NM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("CNRM mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
n2 <- ggplot(subset(NM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
n3 <- ggplot(subset(NM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 
  
g1 <- ggplot(subset(GM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
g2 <- ggplot(subset(GM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
g3 <- ggplot(subset(GM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 

i1 <- ggplot(subset(IM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("IPSL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
i2 <- ggplot(subset(IM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
i3 <- ggplot(subset(IM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 

u1 <- ggplot(subset(UM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("UK mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
u2 <- ggplot(subset(UM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
u3 <- ggplot(subset(UM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 

png(paste0(figp,'corr_ssp585_mesoz_chl_biomes_all_log_LC45.png'), 
    width = 5*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,
           n1,n2,n3,
           g1,g2,g3,
           i1,i2,i3,
           u1,u2,u3,
           nrow = 5, ncol = 3,
           rel_widths = c( 1, 1, 1 ), rel_heights = c( 1, 1, 1 ) ,
           align = 'h' )
dev.off()



### Save corrs as a table
Bcorr <- as.data.frame(t(Ccorr))
Bcorr[3,] <- as.data.frame(t(Ncorr$corr))
Bcorr[4,] <- as.data.frame(t(Gcorr$corr))
Bcorr[5,] <- as.data.frame(t(Icorr$corr))
Bcorr[6,] <- as.data.frame(t(Ucorr$corr))
names(Bcorr) <- c("LC","HCSS","HCPS")
row.names(Bcorr) <- c("biome","CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Bcorr,"SSP585_corrs_mod_chl_biomes_LC45.csv",sep=",",row.names=T)


### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biome==1))

#CNRM
Nlm1 <- lm(Lzmeso~Lchl, data=subset(NM, biome==1))

#GFDL
Glm1 <- lm(Lzmeso~Lchl, data=subset(GM, biome==1))

#IPSL
Ilm1 <- lm(Lzmeso~Lchl, data=subset(IM, biome==1))

#UKESM
Ulm1 <- lm(Lzmeso~Lchl, data=subset(UM, biome==1))

## Save coefficients
Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Nlm1))
Cff1[3,] <- t(coefficients(Glm1))
Cff1[4,] <- t(coefficients(Ilm1))
Cff1[5,] <- t(coefficients(Ulm1))
names(Cff1) <- c("Intercept","Lchl")
row.names(Cff1) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Cff1,"SSP585_coeffs_mod_chl_biome1_LC45.csv",sep=",",row.names=T)

