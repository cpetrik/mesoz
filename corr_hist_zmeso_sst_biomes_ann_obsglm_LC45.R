# Calculate correlation of modeled zmeso biomass and sst
# By model-specific biomes

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Zm <- read.csv("skill_model_obsglm_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv("chl_model_obsglm_all_clim.csv",sep=",",header = T,stringsAsFactors = F)
Sst <- read.csv("sst_model_obsglm_all_clim.csv",sep=",",header = T,stringsAsFactors = F)

CM <- as.data.frame(Zm[,c(1:3,9)])
CM[,5] <- as.data.frame(Chl[,3])
CM[,6] <- as.data.frame(Sst[,3])
names(CM) <- c("Lat","Lon","mesoz","biome","chl","sst")
CM <- na.omit(CM)

NM <- as.data.frame(Zm[,c(1:2,4,10)])
NM[,5] <- as.data.frame(Chl[,4])
NM[,6] <- as.data.frame(Sst[,4])
names(NM) <- c("Lat","Lon","mesoz","biome","chl","sst")
NM <- na.omit(NM)

GM <- as.data.frame(Zm[,c(1:2,5,11)])
GM[,5] <- as.data.frame(Chl[,5])
GM[,6] <- as.data.frame(Sst[,5])
names(GM) <- c("Lat","Lon","mesoz","biome","chl","sst")
GM <- na.omit(GM)

IM <- as.data.frame(Zm[,c(1:2,6,12)])
IM[,5] <- as.data.frame(Chl[,6])
IM[,6] <- as.data.frame(Sst[,6])
names(IM) <- c("Lat","Lon","mesoz","biome","chl","sst")
IM <- na.omit(IM)

UM <- as.data.frame(Zm[,c(1:2,7,13)])
UM[,5] <- as.data.frame(Chl[,7])
UM[,6] <- as.data.frame(Sst[,7])
names(UM) <- c("Lat","Lon","mesoz","biome","chl","sst")
UM <- na.omit(UM)

OG <- as.data.frame(Zm[,c(1:2,8,14)])
OG[,5] <- as.data.frame(Chl[,8])
OG[,6] <- as.data.frame(Sst[,8])
names(OG) <- c("Lat","Lon","mesoz","biome","chl","sst")
OG <- na.omit(OG)

CM$Lzmeso <- log(CM$mesoz+1e-16)
NM$Lzmeso <- log(NM$mesoz+1e-16)
GM$Lzmeso <- log(GM$mesoz+1e-16)
IM$Lzmeso <- log(IM$mesoz+1e-16)
UM$Lzmeso <- log(UM$mesoz+1e-16)
OG$Lzmeso <- log(OG$mesoz+1e-16)


## Restrict latitudes to those observed
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

oid1 <- which(OG$biome==1)
oid2 <- which(OG$Lat > 45 | OG$Lat < -45)
oid <- intersect(oid1,oid2)
OG <- OG[-oid,]

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-sst
Ccorr <- ddply(CM, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))
Ncorr <- ddply(NM, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))
Gcorr <- ddply(GM, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))
Icorr <- ddply(IM, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))
Ucorr <- ddply(UM, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))
Ocorr <- ddply(OG, .(biome), summarize, "corr" = cor(Lzmeso, sst, method = "pearson"))


### Plots

## Zoo & SST biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10()  
c2 <- ggplot(subset(CM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10()  
c3 <- ggplot(subset(CM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10()  

n1 <- ggplot(subset(NM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("CNRM mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10()  
n2 <- ggplot(subset(NM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10()  
n3 <- ggplot(subset(NM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10()  
  
g1 <- ggplot(subset(GM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10()  
g2 <- ggplot(subset(GM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10()  
g3 <- ggplot(subset(GM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10()  

i1 <- ggplot(subset(IM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("IPSL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10()  
i2 <- ggplot(subset(IM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10()  
i3 <- ggplot(subset(IM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10()  

u1 <- ggplot(subset(UM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("UK mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10()  
u2 <- ggplot(subset(UM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10()  
u3 <- ggplot(subset(UM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10()  

o1 <- ggplot(subset(NM,biome==1), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("obs mesoz") + xlab("sst") + #ggtitle("LC") +
  scale_y_log10()  
o2 <- ggplot(subset(NM,biome==2), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("sst") + #ggtitle("HCSS") +
  scale_y_log10()  
o3 <- ggplot(subset(NM,biome==3), aes(y=mesoz, x=sst)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("sst") + #ggtitle("HCPS") +
  scale_y_log10()  

png(paste0(figp,'corr_hist_mesoz_sst_biomes_all_log_obsglm_LC45.png'), 
    width = 5*300,        # 5 x 300 pixels
    height = 9*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,
           n1,n2,n3,
           g1,g2,g3,
           i1,i2,i3,
           u1,u2,u3,
           o1,o2,o3,
           nrow = 6, ncol = 3,
           rel_widths = c( 1, 1, 1 ), rel_heights = c( 1, 1, 1 ) ,
           align = 'h' )
dev.off()



### Save corrs as a table
Bcorr <- as.data.frame(t(Ccorr))
Bcorr[3,] <- as.data.frame(t(Ncorr$corr))
Bcorr[4,] <- as.data.frame(t(Gcorr$corr))
Bcorr[5,] <- as.data.frame(t(Icorr$corr))
Bcorr[6,] <- as.data.frame(t(Ucorr$corr))
Bcorr[7,] <- as.data.frame(t(Ocorr$corr))
names(Bcorr) <- c("LC","HCSS","HCPS")
row.names(Bcorr) <- c("biome","CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Bcorr,"Hist_corrs_mod_sst_biomes_obsglm_LC45.csv",sep=",",row.names=T)


### LM of SST ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm1 <- lm(Lzmeso~sst, data=subset(CM, biome==1))

#CNRM
Nlm1 <- lm(Lzmeso~sst, data=subset(NM, biome==1))

#GFDL
Glm1 <- lm(Lzmeso~sst, data=subset(GM, biome==1))

#IPSL
Ilm1 <- lm(Lzmeso~sst, data=subset(IM, biome==1))

#UKESM
Ulm1 <- lm(Lzmeso~sst, data=subset(UM, biome==1))

#Obs
Olm1 <- lm(Lzmeso~sst, data=subset(OG, biome==1))

## Save coefficients
Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Nlm1))
Cff1[3,] <- t(coefficients(Glm1))
Cff1[4,] <- t(coefficients(Ilm1))
Cff1[5,] <- t(coefficients(Ulm1))
Cff1[6,] <- t(coefficients(Olm1))
names(Cff1) <- c("Intercept","sst")
row.names(Cff1) <- c("CAN","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff1,"Hist_coeffs_mod_sst_biome1_LC_obsglm_LC45.csv",sep=",",row.names=T)

