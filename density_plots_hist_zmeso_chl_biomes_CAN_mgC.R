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

#Units from g to mg
Zm[,3:9] <- 1e3*Zm[,3:9]

CM <- as.data.frame(Zm[,c(1:3,10)])
CM[,5] <- as.data.frame(Chl[,3])
CM[,6] <- as.data.frame(Sst[,3])
names(CM) <- c("Lat","Lon","mesoz","biome","chl","sst")
CM <- na.omit(CM)


OG <- as.data.frame(Zm[,c(1:2,9,16)])
OG[,5] <- as.data.frame(Chl[,9])
OG[,6] <- as.data.frame(Sst[,9])
names(OG) <- c("Lat","Lon","mesoz","biome","chl","sst")
OG <- na.omit(OG)



## Restrict latitudes to those observed
maxlat <- max(OG$Lat)
minlat <- min(OG$Lat)
CM <- subset(CM, Lat >= minlat & Lat <= maxlat)

## Remove LC >45 lat BUT THIS MAY EXPAND OUTSIDE 45DEG UNDER SSP585!?!
cid1 <- which(CM$biome==1)
cid2 <- which(CM$Lat > 45 | CM$Lat < -45)
cid <- intersect(cid1,cid2)
CM <- CM[-cid,]


### Plots -------------------------------------------------------
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl biomes sep
c1 <- ggplot(subset(CM,biome==1), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) +  
  ylab("CAN zmeso") + xlab("chl") + ggtitle("LC") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c2 <- ggplot(subset(CM,biome==2), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + ggtitle("HCSS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c3 <- ggplot(subset(CM,biome==3), aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + ggtitle("HCPS") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')
c0 <- ggplot(CM, aes(y=mesoz, x=chl)) + theme_bw(base_size=12) +  
  geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + 
  geom_smooth(method="lm", se=FALSE, col="red", size = 0.25) + 
  ylab("") + xlab("chl") + ggtitle("Global") +
  scale_y_log10() + scale_x_log10() + theme(legend.position='none')


png(paste0(figp,'scatter_hist_mesoz_chl_biomes_CAN_log_obsglm100_LC45_mgC.png'), 
    width = 6*300,        # 5 x 300 pixels
    height = 2*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,
           nrow = 1, ncol = 3,
           rel_widths = c( 1, 1, 1 ), rel_heights = c( 1, 1, 1 ) ,
           align = 'h' )
dev.off()

png(paste0(figp,'scatter_hist_mesoz_chl_biomes_global_CAN_log_obsglm100_LC45_mgC.png'), 
    width = 8*300,        # 5 x 300 pixels
    height = 2*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c0,
           nrow = 1, ncol = 4,
           rel_widths = c( 1, 1, 1, 1), rel_heights = c( 1, 1, 1, 1 ) ,
           align = 'h' )
dev.off()


