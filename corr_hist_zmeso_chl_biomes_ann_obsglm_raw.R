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

### Data exploration ------------
## Correlations - Jessica said log-log for zoo-chl
Ccorr <- ddply(CM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Mcorr <- ddply(MM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ncorr <- ddply(NM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Gcorr <- ddply(GM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Icorr <- ddply(IM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ucorr <- ddply(UM, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ocorr <- ddply(OG, .(biome), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))


### Plots -------------------------------------------------------
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +  
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

m1 <- ggplot(subset(MM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("CMCC mesoz") + xlab("") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10(limits=c(0.02,0.15)) 
m2 <- ggplot(subset(MM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
m3 <- ggplot(subset(MM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
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

o1 <- ggplot(subset(OM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("obs mesoz") + xlab("chl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() 
o2 <- ggplot(subset(OM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("chl") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() 
o3 <- ggplot(subset(OM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("") + xlab("chl") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() 

png(paste0(figp,'corr_hist_mesoz_chl_biomes_all_log_obsglm_LC45.png'), 
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

png(paste0(figp,'corr_hist_mesoz_chl_biomes_all_log_obsglm_LC45_v2.png'), 
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

png(paste0(figp,'corr_hist_mesoz_chl_biomes_CAN_log_obsglm_LC45.png'), 
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

### Save corrs as a table
Bcorr <- as.data.frame(t(Ccorr))
Bcorr[3,] <- as.data.frame(t(Mcorr$corr))
Bcorr[4,] <- as.data.frame(t(Ncorr$corr))
Bcorr[5,] <- as.data.frame(t(Gcorr$corr))
Bcorr[6,] <- as.data.frame(t(Icorr$corr))
Bcorr[7,] <- as.data.frame(t(Ucorr$corr))
Bcorr[8,] <- as.data.frame(t(Ocorr$corr))
names(Bcorr) <- c("LC","HCSS","HCPS")
row.names(Bcorr) <- c("biome","CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Bcorr,"Hist_corrs_mod_chl_biomes_obsglm_LC45.csv",sep=",",row.names=T)


### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm0 <- lm(Lzmeso~Lchl, data=CM)
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biome==1))
Clm2 <- lm(Lzmeso~Lchl, data=subset(CM, biome==2))
Clm3 <- lm(Lzmeso~Lchl, data=subset(CM, biome==3))
summary(Clm1)
summary(Clm2)
summary(Clm3)

#CMCC
Mlm0 <- lm(Lzmeso~Lchl, data=MM)
Mlm1 <- lm(Lzmeso~Lchl, data=subset(MM, biome==1))
Mlm2 <- lm(Lzmeso~Lchl, data=subset(MM, biome==2))
Mlm3 <- lm(Lzmeso~Lchl, data=subset(MM, biome==3))
summary(Mlm1)
summary(Mlm2)
summary(Mlm3)

#CNRM
Nlm0 <- lm(Lzmeso~Lchl, data=NM)
Nlm1 <- lm(Lzmeso~Lchl, data=subset(NM, biome==1))
Nlm2 <- lm(Lzmeso~Lchl, data=subset(NM, biome==2))
Nlm3 <- lm(Lzmeso~Lchl, data=subset(NM, biome==3))
summary(Nlm1)
summary(Nlm2)
summary(Nlm3)

#GFDL
Glm0 <- lm(Lzmeso~Lchl, data=GM)
Glm1 <- lm(Lzmeso~Lchl, data=subset(GM, biome==1))
Glm2 <- lm(Lzmeso~Lchl, data=subset(GM, biome==2))
Glm3 <- lm(Lzmeso~Lchl, data=subset(GM, biome==3))
summary(Glm1)
summary(Glm2)
summary(Glm3)

#IPSL
Ilm0 <- lm(Lzmeso~Lchl, data=IM)
Ilm1 <- lm(Lzmeso~Lchl, data=subset(IM, biome==1))
Ilm2 <- lm(Lzmeso~Lchl, data=subset(IM, biome==2))
Ilm3 <- lm(Lzmeso~Lchl, data=subset(IM, biome==3))
summary(Ilm1)
summary(Ilm2)
summary(Ilm3)

#UKESM
Ulm0 <- lm(Lzmeso~Lchl, data=UM)
Ulm1 <- lm(Lzmeso~Lchl, data=subset(UM, biome==1))
Ulm2 <- lm(Lzmeso~Lchl, data=subset(UM, biome==2))
Ulm3 <- lm(Lzmeso~Lchl, data=subset(UM, biome==3))
summary(Ulm1)
summary(Ulm2)
summary(Ulm3)

#Obs
Olm0 <- lm(Lzmeso~Lchl, data=OG)
Olm1 <- lm(Lzmeso~Lchl, data=subset(OG, biome==1))
Olm2 <- lm(Lzmeso~Lchl, data=subset(OG, biome==2))
Olm3 <- lm(Lzmeso~Lchl, data=subset(OG, biome==3))
summary(Olm1)
summary(Olm2)
summary(Olm3)

## Save coefficients

Cff0 <- as.data.frame(t(coefficients(Clm0)))
Cff0[2,] <- t(coefficients(Mlm0))
Cff0[3,] <- t(coefficients(Nlm0))
Cff0[4,] <- t(coefficients(Glm0))
Cff0[5,] <- t(coefficients(Ilm0))
Cff0[6,] <- t(coefficients(Ulm0))
Cff0[7,] <- t(coefficients(Olm0))
names(Cff0) <- c("Intercept","Lchl")
row.names(Cff0) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff0,"Hist_coeffs_mod_chl_global_obsglm.csv",sep=",",row.names=T)


Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Mlm1))
Cff1[3,] <- t(coefficients(Nlm1))
Cff1[4,] <- t(coefficients(Glm1))
Cff1[5,] <- t(coefficients(Ilm1))
Cff1[6,] <- t(coefficients(Ulm1))
Cff1[7,] <- t(coefficients(Olm1))
names(Cff1) <- c("Intercept","Lchl")
row.names(Cff1) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff1,"Hist_coeffs_mod_chl_biome1_LC45_obsglm.csv",sep=",",row.names=T)

Cff2 <- as.data.frame(t(coefficients(Clm2)))
Cff2[2,] <- t(coefficients(Mlm2))
Cff2[3,] <- t(coefficients(Nlm2))
Cff2[4,] <- t(coefficients(Glm2))
Cff2[5,] <- t(coefficients(Ilm2))
Cff2[6,] <- t(coefficients(Ulm2))
Cff2[7,] <- t(coefficients(Olm2))
names(Cff2) <- c("Intercept","Lchl")
row.names(Cff2) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff2,"Hist_coeffs_mod_chl_biome2_HCSS_obsglm.csv",sep=",",row.names=T)

Cff3 <- as.data.frame(t(coefficients(Clm3)))
Cff3[2,] <- t(coefficients(Mlm3))
Cff3[3,] <- t(coefficients(Nlm3))
Cff3[4,] <- t(coefficients(Glm3))
Cff3[5,] <- t(coefficients(Ilm3))
Cff3[6,] <- t(coefficients(Ulm3))
Cff3[7,] <- t(coefficients(Olm3))
names(Cff3) <- c("Intercept","Lchl")
row.names(Cff3) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UKESM","obs")
write.table(Cff3,"Hist_coeffs_mod_chl_biome3_HCPS_obsglm.csv",sep=",",row.names=T)

## Save standard error of obsGLMM
se <- as.data.frame(t(sqrt(diag(vcov(Olm0)))))
se[2,] <- sqrt(diag(vcov(Olm1)))
se[3,] <- sqrt(diag(vcov(Olm2)))
se[4,] <- sqrt(diag(vcov(Olm3)))
names(se) <- c("Intercept","Lchl")
row.names(se) <- c("Global","LC","HCSS","HCPS")
write.table(se,"Hist_std_err_biomes_obsglm.csv",sep=",",row.names=T)

## prediction tests
cmin <- log10(0.01)
cmax <- log10(25)
pchl <- exp(log(10)*seq(cmin,cmax,length=50))

plm0 <- predict(Olm0, newdata = data.frame(Lchl = log10(pchl)),
        interval = "confidence")
plm1 <- predict(Olm1, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm2 <- predict(Olm2, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm3 <- predict(Olm3, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")

CI0 <- as.data.frame(plm0)
CI1 <- as.data.frame(plm1)
CI2 <- as.data.frame(plm2)
CI3 <- as.data.frame(plm3)
CI0$Lchl <- log10(pchl)
CI1$Lchl <- log10(pchl)
CI2$Lchl <- log10(pchl)
CI3$Lchl <- log10(pchl)

write.table(CI0,"Hist_CIs_chl_global_obsglm.csv",sep=",",row.names=F)
write.table(CI1,"Hist_CIs_chl_biome1_LC_obsglm.csv",sep=",",row.names=F)
write.table(CI2,"Hist_CIs_chl_biome2_HCSS_obsglm.csv",sep=",",row.names=F)
write.table(CI3,"Hist_CIs_chl_biome3_HCPS_obsglm.csv",sep=",",row.names=F)

