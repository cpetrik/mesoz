# Calculate correlation of modeled zmeso biomass and chl
# By model-specific biomes
# ESMs, obsGLMM100, Stromberg, Moriarty&OBrien
# zmeso was in gC, change to mgC


rm(list=ls())

library(modelsummary)
library(kableExtra)
library(gt)

#source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/")
figp <- "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/"
datap = "/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/"

### load data
## ESM & obsGLMM zmeso, chl
Zm <- read.csv(paste0(datap,"skill_hist_model_obsglm100_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
Chl <- read.csv(paste0(datap,"hist_chl_model_obsglm100_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)
## Stromberg zmeso, chl
SZm <- read.csv(paste0(datap,"skill_hist_Stromberg_all_clim_200.csv"),sep=",",header = T,stringsAsFactors = F)
SChl <- read.csv(paste0(datap,"hist_chl_seawifs_all_clim.csv"),sep=",",header = T,stringsAsFactors = F)
## Moriarty & O'Brien
CZm <- read.csv(paste0(datap,"skill_hist_COPEPOD_all_clim_200_chls.csv"),sep=",",header = T,stringsAsFactors = F)


#Units from g to mg
Zm[,3:9] <- 1e3*Zm[,3:9]
SZm$Stromberg <- 1e3*SZm$Stromberg


#CAN
CM <- as.data.frame(Zm[,c(1:3,10)])
CM[,5] <- as.data.frame(Chl[,3])
names(CM) <- c("Lat","Lon","mesoz","biome","chl")
CM <- na.omit(CM)
#CMCC
MM <- as.data.frame(Zm[,c(1:2,4,11)])
MM[,5] <- as.data.frame(Chl[,4])
names(MM) <- c("Lat","Lon","mesoz","biome","chl")
MM <- na.omit(MM)
#CNRM
NM <- as.data.frame(Zm[,c(1:2,5,12)])
NM[,5] <- as.data.frame(Chl[,5])
names(NM) <- c("Lat","Lon","mesoz","biome","chl")
NM <- na.omit(NM)
#GFDL
GM <- as.data.frame(Zm[,c(1:2,6,13)])
GM[,5] <- as.data.frame(Chl[,6])
names(GM) <- c("Lat","Lon","mesoz","biome","chl")
GM <- na.omit(GM)
#IPSL
IM <- as.data.frame(Zm[,c(1:2,7,14)])
IM[,5] <- as.data.frame(Chl[,7])
names(IM) <- c("Lat","Lon","mesoz","biome","chl")
IM <- na.omit(IM)
#UK
UM <- as.data.frame(Zm[,c(1:2,8,15)])
UM[,5] <- as.data.frame(Chl[,8])
names(UM) <- c("Lat","Lon","mesoz","biome","chl")
UM <- na.omit(UM)
#obsGLMM
OG <- as.data.frame(Zm[,c(1:2,9,16)])
OG[,5] <- as.data.frame(Chl[,9])
names(OG) <- c("Lat","Lon","mesoz","biome","chl")
OG <- na.omit(OG)
#obsSM
OS <- as.data.frame(SZm)
OS[,5] <- as.data.frame(SChl[,3])
names(OS) <- c("Lat","Lon","mesoz","biome","chl")
OS <- na.omit(OS)
#obsMO-S
MOS <- as.data.frame(CZm[,c("Lat","Lon","zmeso200","SEAWIFSchl","SEAWIFSbiomes")])
MOS <- na.omit(MOS)
names(MOS) <- c("Lat","Lon","zmeso200","chl","biome")
#obsMO-M
MOM <- as.data.frame(CZm[,c("Lat","Lon","zmeso200","MODISchl","MODISbiomes")])
MOM <- na.omit(MOM)
names(MOM) <- c("Lat","Lon","zmeso200","chl","biome")


#log10 trans
CM$Lzmeso <- log10(CM$mesoz+1e-16)
MM$Lzmeso <- log10(MM$mesoz+1e-16)
NM$Lzmeso <- log10(NM$mesoz+1e-16)
GM$Lzmeso <- log10(GM$mesoz+1e-16)
IM$Lzmeso <- log10(IM$mesoz+1e-16)
UM$Lzmeso <- log10(UM$mesoz+1e-16)
OG$Lzmeso <- log10(OG$mesoz+1e-16)
OS$Lzmeso <- log10(OS$mesoz+1e-16)
MOS$Lzmeso <- log10(MOS$zmeso200 +1e-16)
MOM$Lzmeso <- log10(MOM$zmeso200 +1e-16)

CM$Lchl <- log10(CM$chl+1e-16)
MM$Lchl <- log10(MM$chl+1e-16)
NM$Lchl <- log10(NM$chl+1e-16)
GM$Lchl <- log10(GM$chl+1e-16)
IM$Lchl <- log10(IM$chl+1e-16)
UM$Lchl <- log10(UM$chl+1e-16)
OG$Lchl <- log10(OG$chl+1e-16)
OS$Lchl <- log10(OS$chl+1e-16)
MOS$Lchl <- log10(MOS$chl+1e-16) 
MOM$Lchl <- log10(MOM$chl+1e-16) 

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

sid1 <- which(OS$biome==1)
sid2 <- which(OS$Lat > 45 | OS$Lat < -45)
sid <- intersect(sid1,sid2)
OS <- OS[-uid,]

sid1 <- which(MOS$biome==1)
sid2 <- which(MOS$Lat > 45 | MOS$Lat < -45)
sid <- intersect(sid1,sid2)
MOS <- MOS[-sid,]

mid1 <- which(MOM$biome==1)
mid2 <- which(MOM$Lat > 45 | MOM$Lat < -45)
mid <- intersect(mid1,mid2)
MOM <- MOM[-mid,]


### LM of Chl ---------------------------------------------------------
##Kwiatkowshki paper does not force through origin

#CAN
Clm0 <- lm(Lzmeso~Lchl, data=CM)
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biome==1))
Clm2 <- lm(Lzmeso~Lchl, data=subset(CM, biome==2))
Clm3 <- lm(Lzmeso~Lchl, data=subset(CM, biome==3))

#CMCC
Mlm0 <- lm(Lzmeso~Lchl, data=MM)
Mlm1 <- lm(Lzmeso~Lchl, data=subset(MM, biome==1))
Mlm2 <- lm(Lzmeso~Lchl, data=subset(MM, biome==2))
Mlm3 <- lm(Lzmeso~Lchl, data=subset(MM, biome==3))

#CNRM
Nlm0 <- lm(Lzmeso~Lchl, data=NM)
Nlm1 <- lm(Lzmeso~Lchl, data=subset(NM, biome==1))
Nlm2 <- lm(Lzmeso~Lchl, data=subset(NM, biome==2))
Nlm3 <- lm(Lzmeso~Lchl, data=subset(NM, biome==3))

#GFDL
Glm0 <- lm(Lzmeso~Lchl, data=GM)
Glm1 <- lm(Lzmeso~Lchl, data=subset(GM, biome==1))
Glm2 <- lm(Lzmeso~Lchl, data=subset(GM, biome==2))
Glm3 <- lm(Lzmeso~Lchl, data=subset(GM, biome==3))

#IPSL
Ilm0 <- lm(Lzmeso~Lchl, data=IM)
Ilm1 <- lm(Lzmeso~Lchl, data=subset(IM, biome==1))
Ilm2 <- lm(Lzmeso~Lchl, data=subset(IM, biome==2))
Ilm3 <- lm(Lzmeso~Lchl, data=subset(IM, biome==3))

#UK
Ulm0 <- lm(Lzmeso~Lchl, data=UM)
Ulm1 <- lm(Lzmeso~Lchl, data=subset(UM, biome==1))
Ulm2 <- lm(Lzmeso~Lchl, data=subset(UM, biome==2))
Ulm3 <- lm(Lzmeso~Lchl, data=subset(UM, biome==3))

#Obs
Olm0 <- lm(Lzmeso~Lchl, data=OG)
Olm1 <- lm(Lzmeso~Lchl, data=subset(OG, biome==1))
Olm2 <- lm(Lzmeso~Lchl, data=subset(OG, biome==2))
Olm3 <- lm(Lzmeso~Lchl, data=subset(OG, biome==3))

#Stromberg 
Slm0 <- lm(Lzmeso~Lchl, data=OS)
Slm1 <- lm(Lzmeso~Lchl, data=subset(OS, biome==1))
Slm2 <- lm(Lzmeso~Lchl, data=subset(OS, biome==2))
Slm3 <- lm(Lzmeso~Lchl, data=subset(OS, biome==3))

#MO Seawifs
CSlm0 <- lm(Lzmeso~Lchl, data=MOS)
CSlm1 <- lm(Lzmeso~Lchl, data=subset(MOS, biome==1))
CSlm2 <- lm(Lzmeso~Lchl, data=subset(MOS, biome==2))
CSlm3 <- lm(Lzmeso~Lchl, data=subset(MOS, biome==3))

#MO MODIS
CMlm0 <- lm(Lzmeso~Lchl, data=MOM)
CMlm1 <- lm(Lzmeso~Lchl, data=subset(MOM, biome==1))
CMlm2 <- lm(Lzmeso~Lchl, data=subset(MOM, biome==2))
CMlm3 <- lm(Lzmeso~Lchl, data=subset(MOM, biome==3))


## Save coefficients
Cff0 <- as.data.frame(t(coefficients(Clm0)))
Cff0[2,] <- t(coefficients(Mlm0))
Cff0[3,] <- t(coefficients(Nlm0))
Cff0[4,] <- t(coefficients(Glm0))
Cff0[5,] <- t(coefficients(Ilm0))
Cff0[6,] <- t(coefficients(Ulm0))
Cff0[7,] <- t(coefficients(Olm0))
Cff0[8,] <- t(coefficients(Slm0))
Cff0[9,] <- t(coefficients(CSlm0))
Cff0[10,] <- t(coefficients(CMlm0))
names(Cff0) <- c("Intercept","Lchl")
row.names(Cff0) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UK","obsGLMM","obsSM","obsMO-S","obsMO-M")
write.table(Cff0,paste0(datap,"Hist_coeffs_mod_chl_global_obsglm100_strom_cope_mgC.csv"),sep=",",row.names=T)


Cff1 <- as.data.frame(t(coefficients(Clm1)))
Cff1[2,] <- t(coefficients(Mlm1))
Cff1[3,] <- t(coefficients(Nlm1))
Cff1[4,] <- t(coefficients(Glm1))
Cff1[5,] <- t(coefficients(Ilm1))
Cff1[6,] <- t(coefficients(Ulm1))
Cff1[7,] <- t(coefficients(Olm1))
Cff1[8,] <- t(coefficients(Slm1))
Cff1[9,] <- t(coefficients(CSlm1))
Cff1[10,] <- t(coefficients(CMlm1))
names(Cff1) <- c("Intercept","Lchl")
row.names(Cff1) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UK","obsGLMM","obsSM","obsMO-S","obsMO-M")
write.table(Cff1,paste0(datap,"Hist_coeffs_mod_chl_biome1_LC45_obsglm100_strom_cope_mgC.csv"),sep=",",row.names=T)

Cff2 <- as.data.frame(t(coefficients(Clm2)))
Cff2[2,] <- t(coefficients(Mlm2))
Cff2[3,] <- t(coefficients(Nlm2))
Cff2[4,] <- t(coefficients(Glm2))
Cff2[5,] <- t(coefficients(Ilm2))
Cff2[6,] <- t(coefficients(Ulm2))
Cff2[7,] <- t(coefficients(Olm2))
Cff2[8,] <- t(coefficients(Slm2))
Cff2[9,] <- t(coefficients(CSlm2))
Cff2[10,] <- t(coefficients(CMlm2))
names(Cff2) <- c("Intercept","Lchl")
row.names(Cff2) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UK","obsGLMM","obsSM","obsMO-S","obsMO-M")
write.table(Cff2,paste0(datap,"Hist_coeffs_mod_chl_biome2_HCSS_obsglm100_strom_cope_mgC.csv"),sep=",",row.names=T)

Cff3 <- as.data.frame(t(coefficients(Clm3)))
Cff3[2,] <- t(coefficients(Mlm3))
Cff3[3,] <- t(coefficients(Nlm3))
Cff3[4,] <- t(coefficients(Glm3))
Cff3[5,] <- t(coefficients(Ilm3))
Cff3[6,] <- t(coefficients(Ulm3))
Cff3[7,] <- t(coefficients(Olm3))
Cff3[8,] <- t(coefficients(Slm3))
Cff3[9,] <- t(coefficients(CSlm3))
Cff3[10,] <- t(coefficients(CMlm3))
names(Cff3) <- c("Intercept","Lchl")
row.names(Cff3) <- c("CAN","CMCC","CNRM","GFDL","IPSL","UK","obsGLMM","obsSM","obsMO-S","obsMO-M")
write.table(Cff3,paste0(datap,"Hist_coeffs_mod_chl_biome3_HCPS_obsglm100_strom_cope_mgC.csv"),sep=",",row.names=T)


### Save all summary stats
## Diff table for each biome
#GLobal
Gmodels <- list(
  "CAN"     = Clm0,
  "CMCC"    = Mlm0,
  "CNRM"    = Nlm0,
  "GFDL"    = Glm0,
  "IPSL"    = Ilm0,
  "UK"      = Ulm0,
  "obsGLMM" = Olm0,
  "obsSM"   = Slm0,
  "obsMO-S" = CSlm0,
  "obsMO-M" = CMlm0
)
modelsummary(Gmodels, fmt=2) #estimate = "stars")
modelsummary(Gmodels, fmt=2, output = paste0(datap,"regress_global_table_mgC.docx"))

#LC
Lmodels <- list(
  "CAN"     = Clm1,
  "CMCC"    = Mlm1,
  "CNRM"    = Nlm1,
  "GFDL"    = Glm1,
  "IPSL"    = Ilm1,
  "UK"      = Ulm1,
  "obsGLMM" = Olm1,
  "obsSM"   = Slm1,
  "obsMO-S" = CSlm1,
  "obsMO-M" = CMlm1
)
modelsummary(Lmodels, fmt=2, output = paste0(datap,"regress_LC_table_mgC.docx"))

#HCSS
Smodels <- list(
  "CAN"     = Clm2,
  "CMCC"    = Mlm2,
  "CNRM"    = Nlm2,
  "GFDL"    = Glm2,
  "IPSL"    = Ilm2,
  "UK"      = Ulm2,
  "obsGLMM" = Olm2,
  "obsSM"   = Slm2,
  "obsMO-S" = CSlm2,
  "obsMO-M" = CMlm2
)
modelsummary(Smodels, fmt=2, output = paste0(datap,"regress_HCSS_table_mgC.docx"))

#HCPS
Pmodels <- list(
  "CAN"     = Clm3,
  "CMCC"    = Mlm3,
  "CNRM"    = Nlm3,
  "GFDL"    = Glm3,
  "IPSL"    = Ilm3,
  "UK"      = Ulm3,
  "obsGLMM" = Olm3,
  "obsSM"   = Slm3,
  "obsMO-S" = CSlm3,
  "obsMO-M" = CMlm3
)
modelsummary(Pmodels, fmt=2, output = paste0(datap,"regress_HCPS_table_mgC.docx"))


