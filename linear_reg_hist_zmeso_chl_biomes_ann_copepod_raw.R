# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes
# COPEPOD data only
# test with 2 diff chls

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

# load data
Zm <- read.csv(paste0(datap,"skill_hist_COPEPOD_all_clim_200_chls.csv"),sep=",",header = T,stringsAsFactors = F)

SM <- as.data.frame(Zm[,c("Lat","Lon","zmeso200","SEAWIFSchl","SEAWIFSbiomes")])
SM <- na.omit(SM)
names(SM) <- c("Lat","Lon","zmeso200","chl","biome")

MM <- as.data.frame(Zm[,c("Lat","Lon","zmeso200","MODISchl","MODISbiomes")])
MM <- na.omit(MM)
names(MM) <- c("Lat","Lon","zmeso200","chl","biome")

# from mgC to gC
SM$zmeso200 <- SM$zmeso200*1e-3
MM$zmeso200 <- MM$zmeso200*1e-3

SM$Lzmeso <- log10(SM$zmeso200+1e-16)
MM$Lzmeso <- log10(MM$zmeso200+1e-16)

SM$Lchl <- log10(SM$chl+1e-16)
MM$Lchl <- log10(MM$chl+1e-16)


## Remove LC >45 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
sid1 <- which(SM$biome==1)
sid2 <- which(SM$Lat > 45 | SM$Lat < -45)
sid <- intersect(sid1,sid2)
SM <- SM[-sid,]

mid1 <- which(MM$biome==1)
mid2 <- which(MM$Lat > 45 | MM$Lat < -45)
mid <- intersect(mid1,mid2)
MM <- MM[-mid,]


### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#Stromberg 
Slm0 <- lm(Lzmeso~Lchl, data=SM)
Slm1 <- lm(Lzmeso~Lchl, data=subset(SM, biome==1))
Slm2 <- lm(Lzmeso~Lchl, data=subset(SM, biome==2))
Slm3 <- lm(Lzmeso~Lchl, data=subset(SM, biome==3))
summary(Slm1)
summary(Slm2)
summary(Slm3)

#Obs
Mlm0 <- lm(Lzmeso~Lchl, data=MM)
Mlm1 <- lm(Lzmeso~Lchl, data=subset(MM, biome==1))
Mlm2 <- lm(Lzmeso~Lchl, data=subset(MM, biome==2))
Mlm3 <- lm(Lzmeso~Lchl, data=subset(MM, biome==3))
summary(Mlm1)
summary(Mlm2)
summary(Mlm3)

## Save coefficients

Cff <- as.data.frame(t(coefficients(Slm1)))
Cff[2,] <- t(coefficients(Slm2))
Cff[3,] <- t(coefficients(Slm3))
Cff[4,] <- t(coefficients(Slm0))
Cff[1,3:4] <- t(coefficients(Mlm1))
Cff[2,3:4] <- t(coefficients(Mlm2))
Cff[3,3:4] <- t(coefficients(Mlm3))
Cff[4,3:4] <- t(coefficients(Mlm0))
names(Cff) <- c("SchlA","SchlB","MchlA","MchlB")
row.names(Cff) <- c("LC","HCSS","HCPS","global")
write.table(Cff,paste0(datap,"Hist_coeffs_chl_global_biomes_copepod.csv"),sep=",",row.names=T)


## Save standard error of obsGLMM
se <- as.data.frame(t(coefficients(Slm0)))
se[2,] <- t(coefficients(Slm1))
se[3,] <- t(coefficients(Slm2))
se[4,] <- t(coefficients(Slm3))

se[1,c(3:4)] <- as.data.frame(t(sqrt(diag(vcov(Slm0)))))
se[2,c(3:4)] <- sqrt(diag(vcov(Slm1)))
se[3,c(3:4)] <- sqrt(diag(vcov(Slm2)))
se[4,c(3:4)] <- sqrt(diag(vcov(Slm3)))
names(se) <- c("aSM","bSM","aSMse","bSMse")
row.names(se) <- c("Global","LC","HCSS","HCPS")
#write.table(se,"Hist_std_err_biomes_stromberg.csv",sep=",",row.names=T)

## prediction tests
cmin <- log10(0.01)
cmax <- log10(25)
pchl <- exp(log(10)*seq(cmin,cmax,length=50))

plm0 <- predict(Slm0, newdata = data.frame(Lchl = log10(pchl)),
        interval = "confidence")
plm1 <- predict(Slm1, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm2 <- predict(Slm2, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")
plm3 <- predict(Slm3, newdata = data.frame(Lchl = log10(pchl)),
                interval = "confidence")

CI0 <- as.data.frame(plm0)
CI1 <- as.data.frame(plm1)
CI2 <- as.data.frame(plm2)
CI3 <- as.data.frame(plm3)
CI0$Lchl <- log10(pchl)
CI1$Lchl <- log10(pchl)
CI2$Lchl <- log10(pchl)
CI3$Lchl <- log10(pchl)

#write.table(CI0,"Hist_CIs_chl_global_stromberg.csv",sep=",",row.names=F)
#write.table(CI1,"Hist_CIs_chl_biome1_LC_stromberg.csv",sep=",",row.names=F)
#write.table(CI2,"Hist_CIs_chl_biome2_HCSS_stromberg.csv",sep=",",row.names=F)
#write.table(CI3,"Hist_CIs_chl_biome3_HCPS_stromberg.csv",sep=",",row.names=F)


m1 <- ggplot(subset(MM,biome==1), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("zmeso") + xlab("obsGLMM chl") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10(limits=c(0.02,0.15)) + theme(legend.position='none')
m2 <- ggplot(subset(MM,biome==2), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m3 <- ggplot(subset(MM,biome==3), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
m0 <- ggplot(MM, aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("") + xlab("") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

s1 <- ggplot(subset(SM,biome==1), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + #, se=FALSE, col="blue", size = 0.25) +  
  ylab("zmeso") + xlab("Stromberg chl") + #ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s2 <- ggplot(subset(SM,biome==2), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s3 <- ggplot(subset(SM,biome==3), aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
s0 <- ggplot(SM, aes(y=zmeso200, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("") + #ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')

png(paste0(figp,'scatter_mesoz_chl_biomes_global_all_log_copepod_LC45.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m0,
           s1,s2,s3,s0,
           nrow = 2, ncol = 4,
           rel_widths = c( 1, 1, 1, 1 ), 
           rel_heights = c(1,1) ,
           align = 'h' )
dev.off()
