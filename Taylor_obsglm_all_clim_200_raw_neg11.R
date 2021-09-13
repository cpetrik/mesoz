# Calculate different skill metrics for each ESM
# log transformed biomass

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


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"

# load data
Tr <- read.csv("skill_model_obsglm_all_clim_200.csv",sep=",",header = T,stringsAsFactors = F)

#----------------------------- Diff Data Scaling ------------------------------------------------
### Standardization 
#1st take log
Tr[,3:8] <- log(Tr[,3:8] + 1e-16)
Tr2 <- Tr[,3:8]

# -1 to 1
Tr11 <- apply(Tr2, MARGIN = 2, FUN = function(X) -1 + 2*((X - min(X,na.rm=T))/diff(range(X,na.rm=T))))
Tr11 <- as.data.frame(Tr11)

### --------------------------- Taylor diagrams -----------------------------
## V1
library(openair)
source(file = "TaylorDiagramCP.r")
source(file = "checkPrep.r")

### Raw (log)
TrC <- Tr2[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- Tr2[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- Tr2[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- Tr2[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- Tr2[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)
Tr3$model <- as.factor(Tr3$model)

cb <- c("#228833","#33BBEE","#004488","#EE6677","#AA3377")
cb1 <- c("#555555","#228833","#33BBEE","#004488","#EE6677","#AA3377")
cb2 <- c("#228833","#33BBEE","#004488","#EE6677","#AA3377","#555555")


t1 <- TaylorDiagram(subset(Tr3,model!="CAN"), obs = "obsGLM", mod = "esm", 
                    group = "model",cols = cb[2:5])
t1C <- TaylorDiagram(subset(Tr3,model=="CAN"), obs = "obsGLM", mod = "esm", 
                    group = "model",cols = cb[1])

pdf( file = paste0(figp,'Taylor_obsglm_all_clim_200_raw.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t1 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
              cols = cb)
dev.off()

# other version
library(plotrix)
source("taylor.diagram.CP.R")
oldpar<-taylor.diagram.CP(Tr2$obsGLM,Tr2$CAN,col=cb[1],pos.cor=FALSE,
                       pcex=1,cex.axis=1,normalize=FALSE)
taylor.diagram(Tr2$obsGLM,Tr2$CNRM,add=TRUE,col=cb[2],pos.cor=FALSE)
taylor.diagram(Tr2$obsGLM,Tr2$GFDL,add=TRUE,col=cb[3],pos.cor=FALSE)
taylor.diagram(Tr2$obsGLM,Tr2$IPSL,add=TRUE,col=cb[4],pos.cor=FALSE)
taylor.diagram(Tr2$obsGLM,Tr2$UK,add=TRUE,col=cb[5],pos.cor=FALSE)
# get approximate legend position
lpos <- 5*sd(Tr2$obsGLM,na.rm = T)
# add a legend
legend(6,13,legend=c("CAN","CNRM","GFDL","IPSL","UK"),pch=19,
       col=cb)




## -1 to 1
TrC <- Tr11[,c("obsGLM","CAN")]
TrC$model <- "CAN"
names(TrC)[2] <- "esm"
TrN <- Tr11[,c("obsGLM","CNRM")]
TrN$model <- "CNRM"
names(TrN)[2] <- "esm"
TrG <- Tr11[,c("obsGLM","GFDL")]
TrG$model <- "GFDL"
names(TrG)[2] <- "esm"
TrI <- Tr11[,c("obsGLM","IPSL")]
TrI$model <- "IPSL"
names(TrI)[2] <- "esm"
TrU <- Tr11[,c("obsGLM","UK")]
TrU$model <- "UK"
names(TrU)[2] <- "esm"
# With CAN
Tr3 <- rbind(TrC,TrN,TrG,TrI,TrU)

pdf( file = paste0(figp,'Taylor_obsglm_all_clim_200_neg11.pdf')) #, width = unit( 5, 'cm' ), height = unit( 5, 'cm' ) )
t5 <- TaylorDiagram(Tr3, obs = "obsGLM", mod = "esm", group = "model",
                    cols = c("green4", "blue", "purple3", "red3", "maroon4"))
dev.off()



