# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
Zm <- read.csv(paste0(datap,"skill_hist_model_obsglm_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv(paste0(datap,"hist_chl_model_obsglm_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)
Sst <- read.csv(paste0(datap,"hist_sst_model_obsglm_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)

CM <- as.data.frame(Zm[,c(1:3,10)])
CM[,5] <- as.data.frame(Chl[,3])
CM[,6] <- as.data.frame(Sst[,3])
names(CM) <- c("Lat","Lon","mesoz","biome","chl","sst")
CM <- na.omit(CM)

MM <- as.data.frame(Zm[,c(1:2,4,11)])
MM[,5] <- as.data.frame(Chl[,4])
MM[,6] <- as.data.frame(Sst[,4])
names(MM) <- c("Lat","Lon","mesoz","biome","chl","sst")
MM <- na.omit(MM)

NM <- as.data.frame(Zm[,c(1:2,5,12)])
NM[,5] <- as.data.frame(Chl[,5])
NM[,6] <- as.data.frame(Sst[,5])
names(NM) <- c("Lat","Lon","mesoz","biome","chl","sst")
NM <- na.omit(NM)

GM <- as.data.frame(Zm[,c(1:2,6,13)])
GM[,5] <- as.data.frame(Chl[,6])
GM[,6] <- as.data.frame(Sst[,6])
names(GM) <- c("Lat","Lon","mesoz","biome","chl","sst")
GM <- na.omit(GM)

IM <- as.data.frame(Zm[,c(1:2,7,14)])
IM[,5] <- as.data.frame(Chl[,7])
IM[,6] <- as.data.frame(Sst[,7])
names(IM) <- c("Lat","Lon","mesoz","biome","chl","sst")
IM <- na.omit(IM)

UM <- as.data.frame(Zm[,c(1:2,8,15)])
UM[,5] <- as.data.frame(Chl[,8])
UM[,6] <- as.data.frame(Sst[,8])
names(UM) <- c("Lat","Lon","mesoz","biome","chl","sst")
UM <- na.omit(UM)

OG <- as.data.frame(Zm[,c(1:2,9,16)])
OG[,5] <- as.data.frame(Chl[,9])
OG[,6] <- as.data.frame(Sst[,9])
names(OG) <- c("Lat","Lon","mesoz","biome","chl","sst")
OG <- na.omit(OG)

CM$Lzmeso <- log10(CM$mesoz+1e-16)
MM$Lzmeso <- log10(MM$mesoz+1e-16)
NM$Lzmeso <- log10(NM$mesoz+1e-16)
GM$Lzmeso <- log10(GM$mesoz+1e-16)
IM$Lzmeso <- log10(IM$mesoz+1e-16)
UM$Lzmeso <- log10(UM$mesoz+1e-16)
OG$Lzmeso <- log10(OG$mesoz+1e-16)

CM$Lchl <- log10(CM$chl+1e-16)
MM$Lchl <- log10(MM$chl+1e-16)
NM$Lchl <- log10(NM$chl+1e-16)
GM$Lchl <- log10(GM$chl+1e-16)
IM$Lchl <- log10(IM$chl+1e-16)
UM$Lchl <- log10(UM$chl+1e-16)
OG$Lchl <- log10(OG$chl+1e-16)

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



### Plots -------------------------------------------------------
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CAN mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c2 <- ggplot(subset(CM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c3 <- ggplot(subset(CM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c0 <- ggplot(CM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

m1 <- ggplot(subset(MM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CMCC mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10(limits=c(0.02,0.15)) + theme(legend.position='none')
m2 <- ggplot(subset(MM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m3 <- ggplot(subset(MM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m0 <- ggplot(MM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 500) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

n1 <- ggplot(subset(NM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CNRM mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n2 <- ggplot(subset(NM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n3 <- ggplot(subset(NM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
n0 <- ggplot(NM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
  
g1 <- ggplot(subset(GM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("GFDL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g2 <- ggplot(subset(GM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g3 <- ggplot(subset(GM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
g0 <- ggplot(GM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

i1 <- ggplot(subset(IM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("IPSL mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i2 <- ggplot(subset(IM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i3 <- ggplot(subset(IM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
i0 <- ggplot(IM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

u1 <- ggplot(subset(UM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + #, se=FALSE, col="blue", size = 0.25) +  
  ylab("UK mesoz") + xlab("") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u2 <- ggplot(subset(UM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u3 <- ggplot(subset(UM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
u0 <- ggplot(UM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

o1 <- ggplot(subset(OG,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("obs mesoz") + xlab("chl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o2 <- ggplot(subset(OG,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("chl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o3 <- ggplot(subset(OG,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("chl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
o0 <- ggplot(OG, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("chl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

png(paste0(figp,'scatter_hist_mesoz_chl_biomes_all_log_obsglm_LC45.png'), 
    width = 5*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,
           n1,n2,n3,
           g1,g2,g3,
           i1,i2,i3,
           u1,u2,u3,
           o1,o2,o3,
           nrow = 6, ncol = 3,
           rel_widths = c( 1, 1, 1 ), rel_heights = c( 1, 1, 1, 1, 1, 1 ) ,
           align = 'h' )
dev.off()

png(paste0(figp,'scatter_hist_mesoz_chl_biomes_all_log_obsglm_LC45_v2.png'), 
    width = 5*300,        # 5 x 300 pixels
    height = 9.25*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,
           NULL,NULL,NULL,
           n1,n2,n3,
           NULL,NULL,NULL,
           g1,g2,g3,
           NULL,NULL,NULL,
           i1,i2,i3,
           NULL,NULL,NULL,
           u1,u2,u3,
           NULL,NULL,NULL,
           o1,o2,o3,
           nrow = 11, ncol = 3,
           rel_widths = c( 1, 1, 1 ), 
           rel_heights = c(1,-0.15,1,-0.15,1,-0.15,1,-0.15,1,-0.15,1) ,
           align = 'h' )
dev.off()

png(paste0(figp,'scatter_hist_mesoz_chl_biomes_CAN_log_obsglm_LC45.png'), 
    width = 6*300,        # 5 x 300 pixels
    height = 2*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1, 1 ), rel_heights = c( 1, 1, 1 ) ,
           align = 'h' )
dev.off()

plot5 <- plot_grid(m1, NULL, m2, rel_heights = c(1, -0.1, 1), align = "hv",
                   labels = c("A", "B"), ncol = 1)
plot5

## With Global
png(paste0(figp,'scatter_hist_mesoz_chl_biomes_global_all_log_obsglm_LC45_v2.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 10*300,
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
           NULL,NULL,NULL,NULL,
           o1,o2,o3,o0,
           nrow = 11, ncol = 4,
           rel_widths = c( 1, 1, 1, 1 ), 
           rel_heights = c(1,-0.15,1,-0.15,1,-0.15,1,-0.15,1,-0.15,1) ,
           align = 'h' )
dev.off()

png(paste0(figp,'scatter_hist_mesoz_chl_biomes_global_CAN_log_obsglm_LC45.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 2*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c0,
           nrow = 1, ncol = 4,
           rel_widths = c( 1, 1, 1, 1), rel_heights = c( 1, 1, 1, 1 ) ,
           align = 'h' )
dev.off()


