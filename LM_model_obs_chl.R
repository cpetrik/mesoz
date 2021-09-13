# Calculate correlation of obs and ESM chl
# Plot obs vs ESM chl
# Violin plot of chl ditr

rm(list=ls())

library(sm)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(PerformanceAnalytics)
library(Hmisc) #rcorr
library(plyr)
library(cowplot) #plot_grid
library(RColorBrewer)
library(Metrics)
library("tidyverse")

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")

# load data
Chl <- read.csv("chl_obs_mod_mean_1950_2014_all_lats.csv",sep=",",header = T,stringsAsFactors = F)

Chl <- na.omit(Chl)

### Data exploration -------------------------------------------------------------------
## Correlations
Smydata <- as.matrix(Chl[,3:8])
Scorr <- rcorr(Smydata) #,type="pearson","spearman"
Scorr_all_p <- data.frame(Scorr$P)
Scorr_all_r <- data.frame(Scorr$r)

## Use log axes instead
Smm <- Chl[,3:8]
Smm <- (10^Smm)
Smm <- as.data.frame(Smm)

### Plots
xlmts2 <- c( 1e-5, 1e-2 ) 
ylmts2 <- c( 3e-3, 10 ) 

m1 <- ggplot(Smm, aes(y=CAN, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CAN chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[2,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[2,1],digits = 2)), size=5)

m2 <- ggplot(Smm, aes(y=CNRM, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("CNRM chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[3,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[3,1],digits = 2)), size=5)

m3 <- ggplot(Smm, aes(y=GFDL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("GFDL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[4,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[4,1],digits = 2)), size=5)

m4 <- ggplot(Smm, aes(y=IPSL, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("IPSL chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[5,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[5,1],digits = 2)), size=5)

m5 <- ggplot(Smm, aes(y=UK, x=obs)) + theme_bw(base_size=14) +  
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + geom_smooth(method="lm", se=FALSE, col="blue", size = 0.75) +  
  ylab("UK chl") + xlab("obs g m-3") + scale_y_log10() + 
  scale_x_log10(limits = xlmts2) + 
  annotate( geom = 'text', y = 1e-3, x = 1e-5, hjust = 0, label=paste0("r = ",signif(Scorr_all_r[6,1],digits = 2)), size=5) +
  annotate( geom = 'text', y = 5e-4, x = 1e-5, hjust = 0, label=paste0("p = ",signif(Scorr_all_p[6,1],digits = 2)), size=5)

pdf( file = 'chl_model_obs_scatter_corr_all_lats.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( m1,m2,m3,m4,m5,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

### linear regression 
cmdl <- lm(Chl$CAN ~ Chl$obs)
cpred <- predict.lmc(mdl)
summary(cmdl)

nmdl <- lm(Chl$CNRM ~ Chl$obs )
npred <- predict.lmc(mdl)
summary(nmdl) 

gmdl <- lm(Chl$GFDL ~ Chl$obs)
gpred <- predict.lmc(mdl)
summary(gmdl) 

imdl <- lm(Chl$IPSL ~ Chl$obs)
ipred <- predict.lmc(mdl)
summary(imdl) 

umdl <- lm(Chl$UK ~ Chl$obs)
upred <- predict.lmc(mdl)
summary(umdl) 

Cff <- as.data.frame(cmdl$coefficients)
Cff[,2] <- nmdl$coefficients
Cff[,3] <- gmdl$coefficients
Cff[,4] <- imdl$coefficients
Cff[,5] <- umdl$coefficients
names(Cff) <- c("CAN","CNRM","GFDL","IPSL","UK")

##
Hdf <- Chl %>%
  select(obs,CAN,CNRM,GFDL,IPSL,UK) %>%
  gather(key = "model", value = "value", -obs)

b1 <- ggplot(Hdf, aes(x = obs, y = value)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + 
  #geom_point(aes(color = model)) + 
  geom_smooth(aes(color = model), method="lm", se=FALSE, size = 0.75) + 
  ggtitle("Hist") + 
  ylab("log10 model chl (g m-3)") + xlab("log10 obs chl (g m-3)") + 
  scale_color_manual(values = c("green4","blue","purple","red","maroon"))


pdf(file = 'chl_model_obs_LM_corr_all_lats.pdf', width = unit( 10, 'cm' ), height = unit( 10, 'cm' ))
b1
dev.off()

pdf( file = 'chl_model_obs_scatterLM_corr_all_lats.pdf', width = unit( 10, 'cm' ), height = unit( 12, 'cm' ) )
plot_grid( m1,m2,m3,m4,m5,b1,
           nrow = 3, ncol = 2,
           rel_widths = c( 1, 1 ), rel_heights = c( 1, 1 ) ,
           align = 'h' )
dev.off()

write.table(Cff,"coeffs_mod_obs_log10chl_all_lat.csv",sep=",",row.names=T)



# ================== Violin plots =======================================
Cdf <- Smm %>%
  select(obs,CAN,CNRM,GFDL,IPSL,UK) %>%
  gather(key = "model", value = "value")

library(scales)
p41 <- ggplot(Cdf, aes(x=model, y=value, fill=model)) + 
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values = c("green4","blue","purple","red","grey","maroon")) +
  stat_summary(fun.y=mean, geom="point", size=2) + 
  labs(title="",x="",y="chl (g m-3)") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() + 
  theme(axis.text.x=element_text(size=12))
p42 <- ggplot(Cdf, aes(x=model, y=value, fill=model)) + 
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values = c("green4","blue","purple","red","grey","maroon")) +
  stat_summary(fun.y=mean, geom="point", size=2) + 
  labs(title="",x="",y="chl (g m-3)") + 
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text.x=element_text(size=12))

p43 <- ggplot(Cdf, aes(x=model, y=value, fill=model)) + 
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values = c("green4","blue","purple","red","grey","maroon")) +
  stat_summary(fun.y=mean, geom="point", size=2) + 
  labs(title="",x="",y="chl (g m-3)") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks()+ theme_bw() + 
  theme(axis.text.x=element_text(size=12))



### 
pdf("model_obs_chl_violin_all_lats.pdf")
p41
dev.off()


