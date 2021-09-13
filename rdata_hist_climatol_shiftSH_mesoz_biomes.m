% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON','units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

czmo_all = zmo_all;
czmo_DJF = zmo_DJF;
czmo_JJA = zmo_JJA;
czmo_MAM = zmo_MAM;
czmo_SON = zmo_SON;

cunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';
load([mpath 'cmcc_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON','units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

mzmo_all = zmo_all;
mzmo_DJF = zmo_DJF;
mzmo_JJA = zmo_JJA;
mzmo_MAM = zmo_MAM;
mzmo_SON = zmo_SON;

munits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON','units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

nzmo_all = zmo_all;
nzmo_DJF = zmo_DJF;
nzmo_JJA = zmo_JJA;
nzmo_MAM = zmo_MAM;
nzmo_SON = zmo_SON;

nunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON')%,'units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

uzmo_all = zmo_all;
uzmo_DJF = zmo_DJF;
uzmo_JJA = zmo_JJA;
uzmo_MAM = zmo_MAM;
uzmo_SON = zmo_SON;

%uunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON');%,'units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

izmo_all = zmo_all;
izmo_DJF = zmo_DJF;
izmo_JJA = zmo_JJA;
izmo_MAM = zmo_MAM;
izmo_SON = zmo_SON;

%iunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat']);

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

gzmo_all = zmo_all;
gzmo_DJF = zmo_DJF;
gzmo_JJA = zmo_JJA;
gzmo_MAM = zmo_MAM;
gzmo_SON = zmo_SON;

%gunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_mesoz.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj] = size(gzmo_all);
glmz_mo = reshape(glmobsmesoz,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);

%% climatologies
ozmo_DJF = nanmean(glmz_mo(:,:,[1 2 12]),3);
ozmo_MAM = nanmean(glmz_mo(:,:,3:5),3);
ozmo_JJA = nanmean(glmz_mo(:,:,6:8),3);
ozmo_SON = nanmean(glmz_mo(:,:,9:11),3);
ozmo_all = nanmean(glmz_mo,3);

%% Biomes
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';

load([fpath 'CanESM_historical_biomes.mat'],'hcbiomes');
load([mpath 'CMCC_historical_biomes.mat'],'hmbiomes');
load([fpath 'CNRM_historical_biomes.mat'],'hnbiomes');
load([fpath 'GFDL_historical_biomes.mat'],'hgbiomes');
load([fpath 'IPSL_historical_biomes.mat'],'hibiomes');
load([fpath 'UKESM_historical_biomes.mat'],'hubiomes','lat','lon');

[lat_b,lon_b] = meshgrid(lat,lon);
clear lat lon

load([bpath 'data_biomes_MODISAqua_x1.mat']);
[lat_m,lon_m] = meshgrid(lat,lon);

%% check orientations
figure(1)
pcolor(czmo_all); shading flat;
title('CAN')

figure(2)
pcolor(mzmo_all); shading flat;
title('CMCC')

figure(7)
pcolor(nzmo_all); shading flat;
title('CNRM')

figure(3)
pcolor(izmo_all); shading flat;
title('IPSL')

figure(4)
pcolor(gzmo_all); shading flat;
title('GFDL')

figure(5)
pcolor(uzmo_all); shading flat;
title('UK')

figure(6)
pcolor(ozmo_all); shading flat;
title('Obs')

%% check biomes
close all
figure(1)
pcolor(hcbiomes); shading flat;
title('CAN')

figure(2)
pcolor(hmbiomes); shading flat;
title('CMCC')

figure(7)
pcolor(hnbiomes); shading flat;
title('CNRM')

figure(3)
pcolor(hibiomes); shading flat;
title('IPSL')

figure(4)
pcolor(hgbiomes); shading flat;
title('GFDL')

figure(5)
pcolor(hubiomes); shading flat;
title('UK')

figure(6)
pcolor(biomes); shading flat;
title('obs')

%% Flip
%UK, CNRM, CAN, CMCC
czmo_all = fliplr(czmo_all);
czmo_DJF = fliplr(czmo_DJF);
czmo_JJA = fliplr(czmo_JJA);
czmo_MAM = fliplr(czmo_MAM);
czmo_SON = fliplr(czmo_SON);

mzmo_all = fliplr(mzmo_all);
mzmo_DJF = fliplr(mzmo_DJF);
mzmo_JJA = fliplr(mzmo_JJA);
mzmo_MAM = fliplr(mzmo_MAM);
mzmo_SON = fliplr(mzmo_SON);

nzmo_all = fliplr(nzmo_all);
nzmo_DJF = fliplr(nzmo_DJF);
nzmo_JJA = fliplr(nzmo_JJA);
nzmo_MAM = fliplr(nzmo_MAM);
nzmo_SON = fliplr(nzmo_SON);

izmo_all = fliplr(izmo_all);
izmo_DJF = fliplr(izmo_DJF);
izmo_JJA = fliplr(izmo_JJA);
izmo_MAM = fliplr(izmo_MAM);
izmo_SON = fliplr(izmo_SON);

uzmo_all = fliplr(uzmo_all);
uzmo_DJF = fliplr(uzmo_DJF);
uzmo_JJA = fliplr(uzmo_JJA);
uzmo_MAM = fliplr(uzmo_MAM);
uzmo_SON = fliplr(uzmo_SON);

biomes = fliplr(biomes);

%% Vectorize, put in ABC order
% convert to same units
%models mol C/m2 -> g C/m2
%obsglm mg C/m2 -> g C/m2
mod_all(:,1) = (czmo_all(:)) * 12.01;
mod_all(:,2) = (mzmo_all(:)) * 12.01;
mod_all(:,3) = (nzmo_all(:)) * 12.01;
mod_all(:,4) = (gzmo_all(:)) * 12.01;
mod_all(:,5) = (izmo_all(:)) * 12.01;
mod_all(:,6) = (uzmo_all(:)) * 12.01;
mod_all(:,7) = (ozmo_all(:)) * 1e-3;

mod_DJF(:,1) = (czmo_DJF(:)) * 12.01;
mod_DJF(:,2) = (mzmo_DJF(:)) * 12.01;
mod_DJF(:,3) = (nzmo_DJF(:)) * 12.01;
mod_DJF(:,4) = (gzmo_DJF(:)) * 12.01;
mod_DJF(:,5) = (izmo_DJF(:)) * 12.01;
mod_DJF(:,6) = (uzmo_DJF(:)) * 12.01;
mod_DJF(:,7) = (ozmo_DJF(:)) * 1e-3;

mod_JJA(:,1) = (czmo_JJA(:)) * 12.01;
mod_JJA(:,2) = (mzmo_JJA(:)) * 12.01;
mod_JJA(:,3) = (nzmo_JJA(:)) * 12.01;
mod_JJA(:,4) = (gzmo_JJA(:)) * 12.01;
mod_JJA(:,5) = (izmo_JJA(:)) * 12.01;
mod_JJA(:,6) = (uzmo_JJA(:)) * 12.01;
mod_JJA(:,7) = (ozmo_JJA(:)) * 1e-3;

mod_MAM(:,1) = (czmo_MAM(:)) * 12.01;
mod_MAM(:,2) = (mzmo_MAM(:)) * 12.01;
mod_MAM(:,3) = (nzmo_MAM(:)) * 12.01;
mod_MAM(:,4) = (gzmo_MAM(:)) * 12.01;
mod_MAM(:,5) = (izmo_MAM(:)) * 12.01;
mod_MAM(:,6) = (uzmo_MAM(:)) * 12.01;
mod_MAM(:,7) = (ozmo_MAM(:)) * 1e-3;

mod_SON(:,1) = (czmo_SON(:)) * 12.01;
mod_SON(:,2) = (mzmo_SON(:)) * 12.01;
mod_SON(:,3) = (nzmo_SON(:)) * 12.01;
mod_SON(:,4) = (gzmo_SON(:)) * 12.01;
mod_SON(:,5) = (izmo_SON(:)) * 12.01;
mod_SON(:,6) = (uzmo_SON(:)) * 12.01;
mod_SON(:,7) = (ozmo_SON(:)) * 1e-3;

%% ADD BIOMES
bvec(:,1) = hcbiomes(:);
bvec(:,2) = hmbiomes(:);
bvec(:,3) = hnbiomes(:);
bvec(:,4) = hgbiomes(:);
bvec(:,5) = hibiomes(:);
bvec(:,6) = hubiomes(:);
bvec(:,7) = biomes(:);

%% all clim
lat = lat_g(:);
lon = lon_g(:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:9) = mod_all;
comb(:,10:16) = bvec;

% nn = ~isnan(comb(:,8));
% comb = comb(nn,:);

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome','obsbiome'});
writetable(obsmod,'skill_hist_model_obsglm_all_clim_200.csv')

%% Winter
lat = lat_g(:);
lon = lon_g(:);
dcomb(:,1) = lat;
dcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
dcomb(nh,3:9) = mod_DJF(nh,:);
dcomb(sh,3:9) = mod_JJA(sh,:);
dcomb(:,10:16) = bvec;

nn = ~isnan(dcomb(:,8));
dcomb = dcomb(nn,:);

omDJF = array2table(dcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome','obsbiome'});
writetable(omDJF,'skill_hist_model_obsglm_DJF_clim_200.csv')

%% Summer clim
lat = lat_g(:);
lon = lon_g(:);

jcomb(:,1) = lat;
jcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
jcomb(nh,3:9) = mod_JJA(nh,:);
jcomb(sh,3:9) = mod_DJF(sh,:);
jcomb(:,10:16) = bvec;

nn = ~isnan(jcomb(:,8));
jcomb = jcomb(nn,:);

omJJA = array2table(jcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome','obsbiome'});
writetable(omJJA,'skill_hist_model_obsglm_JJA_clim_200.csv')

%% Spring clim
lat = lat_g(:);
lon = lon_g(:);

mcomb(:,1) = lat;
mcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
mcomb(nh,3:9) = mod_MAM(nh,:);
mcomb(sh,3:9) = mod_SON(sh,:);
mcomb(:,10:16) = bvec;

nn = ~isnan(mcomb(:,8));
mcomb = mcomb(nn,:);

omMAM = array2table(mcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome','obsbiome'});
writetable(omMAM,'skill_hist_model_obsglm_MAM_clim_200.csv')

%% Fall clim
lat = lat_g(:);
lon = lon_g(:);

scomb(:,1) = lat;
scomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
scomb(nh,3:9) = mod_SON(nh,:);
scomb(sh,3:9) = mod_MAM(sh,:);
scomb(:,10:16) = bvec;

nn = ~isnan(scomb(:,8));
scomb = scomb(nn,:);

omSON = array2table(scomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome','obsbiome'});
writetable(omSON,'skill_hist_model_obsglm_SON_clim_200.csv')

%%
save('skill_hist_model_obsglm_climatols.mat','comb','dcomb','jcomb','mcomb','scomb')

