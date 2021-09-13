% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';
load([cpath 'can_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],...
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
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/ssp585/';
load([mpath 'cmcc_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],...
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
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';
load([npath 'cnrm_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],...
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
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
load([upath 'ukesm_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],...
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
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
load([ipath 'ipsl_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],...
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
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
load([gpath 'gfdl_ssp585_zmeso200_onedeg_climatol_2051_2100.mat']);

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

[lat_g,lon_g] = meshgrid(lat,lon);

%gunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% Biomes
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';

load([fpath 'CanESM_ssp585_biomes.mat'],'scbiomes');
load([mpath 'CMCC_ssp585_biomes.mat'],'smbiomes');
load([fpath 'CNRM_ssp585_biomes.mat'],'snbiomes');
load([fpath 'GFDL_ssp585_biomes.mat'],'sgbiomes');
load([fpath 'IPSL_ssp585_biomes.mat'],'sibiomes');
load([fpath 'UKESM_ssp585_biomes.mat'],'subiomes','lat','lon');

[lat_b,lon_b] = meshgrid(lat,lon);
clear lat lon

%% check orientations
figure(1)
pcolor(czmo_all); shading flat;
title('CAN')

figure(7)
pcolor(mzmo_all); shading flat;
title('CMCC')

figure(2)
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
pcolor(lat_g); shading flat;
title('latg')

%% check orientations
close all
figure(1)
pcolor(scbiomes); shading flat;
title('CAN')

figure(7)
pcolor(smbiomes); shading flat;
title('CMCC')

figure(2)
pcolor(snbiomes); shading flat;
title('CNRM')

figure(3)
pcolor(sibiomes); shading flat;
title('IPSL')

figure(4)
pcolor(sgbiomes); shading flat;
title('GFDL')

figure(5)
pcolor(subiomes); shading flat;
title('UK')

figure(6)
pcolor(lat_b); shading flat;
title('latb')


%% Flip
%CMCC, CNRM, CAN, IPSL
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

mod_DJF(:,1) = (czmo_DJF(:)) * 12.01;
mod_DJF(:,2) = (mzmo_DJF(:)) * 12.01;
mod_DJF(:,3) = (nzmo_DJF(:)) * 12.01;
mod_DJF(:,4) = (gzmo_DJF(:)) * 12.01;
mod_DJF(:,5) = (izmo_DJF(:)) * 12.01;
mod_DJF(:,6) = (uzmo_DJF(:)) * 12.01;

mod_JJA(:,1) = (czmo_JJA(:)) * 12.01;
mod_JJA(:,2) = (mzmo_JJA(:)) * 12.01;
mod_JJA(:,3) = (nzmo_JJA(:)) * 12.01;
mod_JJA(:,4) = (gzmo_JJA(:)) * 12.01;
mod_JJA(:,5) = (izmo_JJA(:)) * 12.01;
mod_JJA(:,6) = (uzmo_JJA(:)) * 12.01;

mod_MAM(:,1) = (czmo_MAM(:)) * 12.01;
mod_MAM(:,2) = (mzmo_MAM(:)) * 12.01;
mod_MAM(:,3) = (nzmo_MAM(:)) * 12.01;
mod_MAM(:,4) = (gzmo_MAM(:)) * 12.01;
mod_MAM(:,5) = (izmo_MAM(:)) * 12.01;
mod_MAM(:,6) = (uzmo_MAM(:)) * 12.01;

mod_SON(:,1) = (czmo_SON(:)) * 12.01;
mod_SON(:,2) = (mzmo_SON(:)) * 12.01;
mod_SON(:,3) = (nzmo_SON(:)) * 12.01;
mod_SON(:,4) = (gzmo_SON(:)) * 12.01;
mod_SON(:,5) = (izmo_SON(:)) * 12.01;
mod_SON(:,6) = (uzmo_SON(:)) * 12.01;

%% ADD BIOMES
bvec(:,1) = scbiomes(:);
bvec(:,2) = smbiomes(:);
bvec(:,3) = snbiomes(:);
bvec(:,4) = sgbiomes(:);
bvec(:,5) = sibiomes(:);
bvec(:,6) = subiomes(:);

%% all clim
lat = lat_g(:);
lon = lon_g(:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:8) = mod_all;
comb(:,9:14) = bvec;

% nn = ~isnan(comb(:,1)); % Need matching mesoz-chl per model for corrs
% comb = comb(nn,:);

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome'});
writetable(obsmod,'ssp585_mesoz_all_clim_200.csv')

%% Winter
lat = lat_g(:);
lon = lon_g(:);
dcomb(:,1) = lat;
dcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
dcomb(nh,3:8) = mod_DJF(nh,:);
dcomb(sh,3:8) = mod_JJA(sh,:);
dcomb(:,9:14) = bvec;

omDJF = array2table(dcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome'});
writetable(omDJF,'ssp585_mesoz_DJF_clim_200.csv')

%% Summer clim
lat = lat_g(:);
lon = lon_g(:);

jcomb(:,1) = lat;
jcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
jcomb(nh,3:8) = mod_JJA(nh,:);
jcomb(sh,3:8) = mod_DJF(sh,:);
jcomb(:,9:14) = bvec;

omJJA = array2table(jcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome'});
writetable(omJJA,'ssp585_mesoz_JJA_clim_200.csv')

%% Spring clim
lat = lat_g(:);
lon = lon_g(:);

mcomb(:,1) = lat;
mcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
mcomb(nh,3:8) = mod_MAM(nh,:);
mcomb(sh,3:8) = mod_SON(sh,:);
mcomb(:,9:14) = bvec;

omMAM = array2table(mcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome'});
writetable(omMAM,'ssp585_mesoz_MAM_clim_200.csv')

%% Fall clim
lat = lat_g(:);
lon = lon_g(:);

scomb(:,1) = lat;
scomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
scomb(nh,3:8) = mod_SON(nh,:);
scomb(sh,3:8) = mod_MAM(sh,:);
scomb(:,9:14) = bvec;

omSON = array2table(scomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK',...
    'CAbiome','CMbiome','CNbiome','GFbiome','IPbiome','UKbiome'});
writetable(omSON,'ssp585_mesoz_SON_clim_200.csv')

%%
save('ssp585_mesoz_climatols.mat','comb','dcomb','jcomb','mcomb','scomb')
