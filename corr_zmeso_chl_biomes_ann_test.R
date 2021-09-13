# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes

rm(list=ls())

library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
CM <- read.csv("can_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)
### GET OBS CHL and SST w/BIOMES !
#Ann <- read.csv("obs_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv",sep=",",header = T,stringsAsFactors = F)

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
CM$chl <- CM$chl * 1e3
GM$chl <- GM$chl * 1e3
UM$chl <- UM$chl * 1e3

CM$Lzmeso <- log(CM$zmeso+1e-16)
NM$Lzmeso <- log(NM$zmeso+1e-16)
GM$Lzmeso <- log(GM$zmeso+1e-16)
IM$Lzmeso <- log(IM$zmeso+1e-16)
UM$Lzmeso <- log(UM$zmeso+1e-16)

CM$Lchl <- log(CM$chl+1e-16)
NM$Lchl <- log(NM$chl+1e-16)
GM$Lchl <- log(GM$chl+1e-16)
IM$Lchl <- log(IM$chl+1e-16)
UM$Lchl <- log(UM$chl+1e-16)

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl

Ccorr <- ddply(CM, .(biomes), summarize, "corr" = cor(log(zmeso+1e-16), log(chl+1e-16), 
              method = "pearson"))
Ccorr2 <- ddply(CM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))

###
library(Hmisc)
corrByGroup <- function(xx){
  return(data.frame(cbind(correl = round(rcorr(xx$Lzmeso, xx$Lchl)$r[1,2], digits=3),
                          n = rcorr(xx$Lzmeso, xx$Lchl)$n[1,2],
                          pvalue = round(rcorr(xx$Lzmeso, xx$Lchl)$P[1,2], digits=3))))
}
Ccorr3 <- corrByGroup(CM) #not divided by biome

###
library(tidyverse)
library(broom)

CM  %>% 
  group_by(biomes) %>%
  summarize(correlation = cor(Lzmeso, Lchl, method = "pearson"))
# A tibble: 3 x 2
# biomes correlation - same result as Ccorr and Ccorr2 methods
# <int>       <dbl>
# 1      1       0.439
# 2      2       0.531
# 3      3       0.625

# with pvalues and further stats
CM %>% 
  nest(-biomes) %>% 
  mutate(cor=map(data,~cor.test(.x$Lzmeso, .x$Lchl, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T)
# A tibble: 3 x 11
# biomes data      cor   estimate statistic   p.value parameter conf.low conf.high
# <int> <list>    <lis>    <dbl>     <dbl>     <dbl>     <int>    <dbl>     <dbl>
# 1      1 <tibble … <hte…    0.439      84.7 0.            30039    0.430     0.448
# 2      3 <tibble … <hte…    0.625      67.3 0.             7059    0.611     0.639
# 3      2 <tibble … <hte…    0.531      39.8 6.10e-293      4043    0.508     0.553

Ncorr <- ddply(NM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Gcorr <- ddply(GM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Icorr <- ddply(IM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
Ucorr <- ddply(UM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))
#Acorr <- ddply(AM, .(biomes), summarize, "corr" = cor(Lzmeso, Lchl, method = "pearson"))


### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl all biomes together
a1 <- ggplot(CM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") + #, se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() #+ 
   # annotate( geom = 'text', y = 1e3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
   # annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(NM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a3 <- ggplot(GM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a4 <- ggplot(IM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
a5 <- ggplot(UM, aes(y=zmeso, x=chl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() 
  
# untransform observations
Amz$Robs <- 10^Amz$obs
Amz$Rchl <- 10^Amz$chl * 1e-3 #from mg to g

a6 <- ggplot(Amz, aes(y=Robs, x=Rchl, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm") +  
  ylab("Obs zmeso mgC m-2") + xlab("Obs chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Acorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_hist_mesoz_chl_biomes_all_log.pdf', width = unit( 12, 'cm' ), height = unit( 10, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,#a6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 1 LC
b1 <- ggplot(subset(CM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr[1,2],digits = 2)), size=5) 

b2 <- ggplot(subset(NM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr[1,2],digits = 2)), size=5) 

b3 <- ggplot(subset(GM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr[1,2],digits = 2)), size=5) 

b4 <- ggplot(subset(IM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr[1,2],digits = 2)), size=5) 

b5 <- ggplot(subset(UM, biomes==1), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e1, x = 3e-6, hjust = 0, label=paste0("r = ",signif(Ucorr[1,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_chl_biome1_LC_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( b1,b2,b3,b4,b5,#b6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 2 HCSS
c1 <- ggplot(subset(CM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[2,2],digits = 2)), size=5) 

c2 <- ggplot(subset(NM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[2,2],digits = 2)), size=5) 

c3 <- ggplot(subset(GM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[2,2],digits = 2)), size=5) 

c4 <- ggplot(subset(IM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[2,2],digits = 2)), size=5) 

c5 <- ggplot(subset(UM, biomes==2), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[2,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_chl_biome2_HCSS_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( c1,c2,c3,c4,c5,#c6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & Chl Biome 3 HCPS
d1 <- ggplot(subset(CM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[3,2],digits = 2)), size=5) 

d2 <- ggplot(subset(NM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[3,2],digits = 2)), size=5) 

d3 <- ggplot(subset(GM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[3,2],digits = 2)), size=5) 

d4 <- ggplot(subset(IM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[3,2],digits = 2)), size=5) 

d5 <- ggplot(subset(UM, biomes==3), aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl g m-3") + 
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[3,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_chl_biome3_HCPS_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( d1,d2,d3,d4,d5,#d6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(CM, aes(y=zmeso, x=sst, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") 

s2 <- ggplot(NM, aes(y=zmeso, x=sst, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") 

s3 <- ggplot(GM, aes(y=zmeso, x=sst, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") 

s4 <- ggplot(IM, aes(y=zmeso, x=sst, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") 

s5 <- ggplot(UM, aes(y=zmeso, x=sst, color=factor(biomes))) + theme_bw(base_size=14) +  
  geom_point(alpha=0.5) + geom_smooth(method="lm") +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") 

pdf( file = 'corr_hist_mesoz_sst_biomes_all_log.pdf', width = unit( 10, 'cm' ), height = unit( 10, 'cm' ) )
plot_grid( s1,s2,s3,s4,s5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST Biome 1 LC
b1 <- ggplot(subset(CM, biomes==1), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e2, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ccorr[1,2],digits = 2)), size=5) 

b2 <- ggplot(subset(NM, biomes==1), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr[1,2],digits = 2)), size=5) 

b3 <- ggplot(subset(GM, biomes==1), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr[1,2],digits = 2)), size=5) 

b4 <- ggplot(subset(IM, biomes==1), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr[1,2],digits = 2)), size=5) 

b5 <- ggplot(subset(UM, biomes==1), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="red", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e1, x = 3e-6, hjust = 0, label=paste0("r = ",signif(Ucorr[1,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_sst_biome1_LC_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( b1,b2,b3,b4,b5,#b6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST Biome 2 HCSS
c1 <- ggplot(subset(CM, biomes==2), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[2,2],digits = 2)), size=5) 

c2 <- ggplot(subset(NM, biomes==2), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[2,2],digits = 2)), size=5) 

c3 <- ggplot(subset(GM, biomes==2), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 2e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[2,2],digits = 2)), size=5) 

c4 <- ggplot(subset(IM, biomes==2), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[2,2],digits = 2)), size=5) 

c5 <- ggplot(subset(UM, biomes==2), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="green4", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[2,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_sst_biome2_HCSS_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( c1,c2,c3,c4,c5,#c6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST Biome 3 HCPS
d1 <- ggplot(subset(CM, biomes==3), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ccorr[3,2],digits = 2)), size=5) 

d2 <- ggplot(subset(NM, biomes==3), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ncorr[3,2],digits = 2)), size=5) 

d3 <- ggplot(subset(GM, biomes==3), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Gcorr[3,2],digits = 2)), size=5) 

d4 <- ggplot(subset(IM, biomes==3), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) + 
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Icorr[3,2],digits = 2)), size=5) 

d5 <- ggplot(subset(UM, biomes==3), aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  scale_y_log10() +  
  annotate( geom = 'text', y = 1e0, x = 1e-4, hjust = 0, label=paste0("r = ",signif(Ucorr[3,2],digits = 2)), size=5) 

pdf( file = 'corr_hist_mesoz_sst_biome3_HCPS_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( d1,d2,d3,d4,d5,#d6,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()




### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==1))
Clm2 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==2))
Clm3 <- lm(Lzmeso~Lchl, data=subset(CM, biomes==3))
summary(Clm1)
summary(Clm2)
summary(Clm3)

#CNRM
Nlm1 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==1))
Nlm2 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==2))
Nlm3 <- lm(Lzmeso~Lchl, data=subset(NM, biomes==3))
summary(Nlm1)
summary(Nlm2)
summary(Nlm3)

#GFDL
Glm1 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==1))
Glm2 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==2))
Glm3 <- lm(Lzmeso~Lchl, data=subset(GM, biomes==3))
summary(Glm1)
summary(Glm2)
summary(Glm3)

#IPSL
Ilm1 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==1))
Ilm2 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==2))
Ilm3 <- lm(Lzmeso~Lchl, data=subset(IM, biomes==3))
summary(Ilm1)
summary(Ilm2)
summary(Ilm3)

#UKESM
Ulm1 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==1))
Ulm2 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==2))
Ulm3 <- lm(Lzmeso~Lchl, data=subset(UM, biomes==3))
summary(Ulm1)
summary(Ulm2)
summary(Ulm3)

#Obs
Olm1 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==1))
Olm2 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==2))
Olm3 <- lm(Lzmeso~Lchl, data=subset(OM, biomes==3))
summary(Olm1)
summary(Olm2)
summary(Olm3)

## Save coefficients
Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Nlm1))
Cff1[3,] <- t(coefficients(Glm1))
Cff1[4,] <- t(coefficients(Ilm1))
Cff1[5,] <- t(coefficients(Ulm1))
#Cff1[6,] <- t(coefficients(Olm1))
names(Cff1) <- c("Intercept","Lchl")
row.names(Cff1) <- c("CAN","CNRM","GFDL","IPSL","UKESM") #,"obs")
write.table(Cff1,"Hist_coeffs_mod_chl_biome1_LC.csv",sep=",",row.names=T)

Cff2 <- as.data.frame(t(coefficients(Clm2)))
Cff2[2,] <- t(coefficients(Nlm2))
Cff2[3,] <- t(coefficients(Glm2))
Cff2[4,] <- t(coefficients(Ilm2))
Cff2[5,] <- t(coefficients(Ulm2))
#Cff2[6,] <- t(coefficients(Olm2))
names(Cff2) <- c("Intercept","Lchl")
row.names(Cff2) <- c("CAN","CNRM","GFDL","IPSL","UKESM") #,"obs")
write.table(Cff2,"Hist_coeffs_mod_chl_biome2_HCSS.csv",sep=",",row.names=T)

Cff3 <- as.data.frame(t(coefficients(Clm3)))
Cff3[2,] <- t(coefficients(Nlm3))
Cff3[3,] <- t(coefficients(Glm3))
Cff3[4,] <- t(coefficients(Ilm3))
Cff3[5,] <- t(coefficients(Ulm3))
#Cff3[6,] <- t(coefficients(Olm3))
names(Cff3) <- c("Intercept","Lchl")
row.names(Cff3) <- c("CAN","CNRM","GFDL","IPSL","UKESM") #,"obs")
write.table(Cff3,"Hist_coeffs_mod_chl_biome3_HCPS.csv",sep=",",row.names=T)



### LM of SST ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin
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

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"Hist_coeffs_mod_sst_trop_all_clim_200.csv",sep=",",row.names=T)


# =========================================================================
#  SSP 585 
# =========================================================================

# load data
CM <- read.csv("can_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
NM <- read.csv("cnrm_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
GM <- read.csv("gfdl_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
IM <- read.csv("ipsl_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)
UM <- read.csv("ukesm_ssp585_mod_zoo_chl_sst_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

### Compare just to growing season
# Latitude >= 60 N: JJA
# Latitude <= -60 S: DJF
# Latitude < 60 N & >= 30 N: JJA + MAM 
# Latitude > - 60 N & <= -30 N: SON + DJF 
# Latitude between 30 N and 30 S: all seasons

## All seasons
aCM <- subset(CM, Lat > -30 & Lat < 30)
aNM <- subset(NM, Lat > -30 & Lat < 30)
aGM <- subset(GM, Lat > -30 & Lat < 30)
aIM <- subset(IM, Lat > -30 & Lat < 30)
aUM <- subset(UM, Lat > -30 & Lat < 30)

### Correct units, set all to g/m3
#CNRM and IPSL accidentally uploaded as g/m3, others all kg/m3
aCM$chl <- aCM$chl * 1e3
aGM$chl <- aGM$chl * 1e3
aUM$chl <- aUM$chl * 1e3

### Data exploration -------------------------------------------------------------------
## Correlations - Jessica said log-log for zoo-chl
Cmydata <- as.matrix(aCM[,3:5])
Cmydata[,1:2] <- log10(Cmydata[,1:2])
Ccorr <- rcorr(Cmydata) #,type="pearson","spearman"
Ccorr_all_p <- data.frame(Ccorr$P)
Ccorr_all_r <- data.frame(Ccorr$r)

Nmydata <- as.matrix(aNM[,3:5])
Nmydata[,1:2] <- log10(Nmydata[,1:2])
Ncorr <- rcorr(Nmydata) #,type="pearson","spearman"
Ncorr_all_p <- data.frame(Ncorr$P)
Ncorr_all_r <- data.frame(Ncorr$r)

Gmydata <- as.matrix(aGM[,3:5])
Gmydata[,1:2] <- log10(Gmydata[,1:2])
Gcorr <- rcorr(Gmydata) #,type="pearson","spearman"
Gcorr_all_p <- data.frame(Gcorr$P)
Gcorr_all_r <- data.frame(Gcorr$r)

Imydata <- as.matrix(aIM[,3:5])
Imydata[,1:2] <- log10(Imydata[,1:2])
Icorr <- rcorr(Imydata) #,type="pearson","spearman"
Icorr_all_p <- data.frame(Icorr$P)
Icorr_all_r <- data.frame(Icorr$r)

Umydata <- as.matrix(aUM[,3:5])
Umydata[,1:2] <- log10(Umydata[,1:2])
Ucorr <- rcorr(Umydata) #,type="pearson","spearman"
Ucorr_all_p <- data.frame(Ucorr$P)
Ucorr_all_r <- data.frame(Ucorr$r)



### Plots
xlmts2 <- c( 1, 3e4 ) 

## Zoo & Chl all biomes together
a1 <- ggplot(aCM, aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN chl kg m-3") + 
  scale_y_log10() + scale_x_log10() + 
  # annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 1e2, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[2,1],digits = 2)), size=5)

a2 <- ggplot(aNM, aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 8e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[2,1],digits = 2)), size=5)

a3 <- ggplot(aGM, aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 3e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 2e3, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[2,1],digits = 2)), size=5)

a4 <- ggplot(aIM, aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL chl kg m-3") + 
  # annotate( geom = 'text', y = 4000, x = 0, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 3500, x = 0, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 1e3, x = 3e-5, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 700, x = 3e-5, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[2,1],digits = 2)), size=5)

a5 <- ggplot(aUM, aes(y=zmeso, x=chl)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM chl kg m-3") + 
  # annotate( geom = 'text', y = 6000, x = 0, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  # annotate( geom = 'text', y = 5500, x = 0, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)
  scale_y_log10() + scale_x_log10() + 
  annotate( geom = 'text', y = 2e4, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 1e4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[2,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_chl_tropics_ann_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( a1,a2,a3,a4,a5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()


## Zoo & SST
s1 <- ggplot(aCM, aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN zmeso mgC m-2") + xlab("CAN sst") + 
  annotate( geom = 'text', y = 600, x = 15, hjust = 0, label=paste0("r = ",signif(Ccorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 550, x = 15, hjust = 0, label=paste0("p = ",signif(Ccorr_all_p[3,1],digits = 2)), size=5)
# scale_y_log10() + scale_x_log10() + 
# annotate( geom = 'text', y = 1e-8, x = 2e3, hjust = 0, label=paste0("r = ",signif(Acorr_all_r[3,1],digits = 2)), size=5) +
# annotate( geom = 'text', y = 1e-9, x = 2e3, hjust = 0, label=paste0("rmse = ",signif(aRMSD[3,1],digits = 2)), size=5)

s2 <- ggplot(aNM, aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM zmeso mgC m-2") + xlab("CNRM sst") + 
  annotate( geom = 'text', y = 6000, x = 20, hjust = 0, label=paste0("r = ",signif(Ncorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 20, hjust = 0, label=paste0("p = ",signif(Ncorr_all_p[3,1],digits = 2)), size=5)

s3 <- ggplot(aGM, aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL zmeso mgC m-2") + xlab("GFDL sst") + 
  annotate( geom = 'text', y = 3300, x = 17, hjust = 0, label=paste0("r = ",signif(Gcorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 3000, x = 17, hjust = 0, label=paste0("p = ",signif(Gcorr_all_p[3,1],digits = 2)), size=5)

s4 <- ggplot(aIM, aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL zmeso mgC m-2") + xlab("IPSL sst") + 
  annotate( geom = 'text', y = 1000, x = 15, hjust = 0, label=paste0("r = ",signif(Icorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 900, x = 15, hjust = 0, label=paste0("p = ",signif(Icorr_all_p[3,1],digits = 2)), size=5)

s5 <- ggplot(aUM, aes(y=zmeso, x=sst)) + theme_bw(base_size=14) +  
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UKESM zmeso mgC m-2") + xlab("UKESM sst") + 
  annotate( geom = 'text', y = 6000, x = 18, hjust = 0, label=paste0("r = ",signif(Ucorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5500, x = 18, hjust = 0, label=paste0("p = ",signif(Ucorr_all_p[3,1],digits = 2)), size=5)

pdf( file = 'corr_ssp585_mesoz_sst_tropics_ann_clim_200_log.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
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
write.table(Cff,"SSP585_coeffs_mod_chl_trop_all_clim_200.csv",sep=",",row.names=T)

Sff <- as.data.frame(t(coefficients(ClmS)))
Sff[2,] <- t(coefficients(NlmS))
Sff[3,] <- t(coefficients(GlmS))
Sff[4,] <- t(coefficients(IlmS))
Sff[5,] <- t(coefficients(UlmS))
row.names(Sff) <- c("CAN","CNRM","GFDL","IPSL","UKESM")
write.table(Sff,"SSP585_coeffs_mod_sst_trop_all_200.csv",sep=",",row.names=T)

