% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% depth info
Cdir = '/Volumes/MIP/Fish-MIP/CMIP6/';
load([Cdir 'GFDL/gridspec_gfdl_cmip6.mat'],'deptho','LAT','LON','lmask');

%% COPEPOD ---------------------------------------
load('copepod-2012_cmass_all_gridded.mat',...
    'zoo_g','fileid','units')
cope_all = zoo_g;
clear zoo_g

load('copepod-2012_cmass_DJF_gridded.mat','zoo_g')
cope_DJF = zoo_g;
clear zoo_g

load('copepod-2012_cmass_MAM_gridded.mat','zoo_g')
cope_MAM = zoo_g;
clear zoo_g

load('copepod-2012_cmass_JJA_gridded.mat','zoo_g')
cope_JJA = zoo_g;
clear zoo_g

load('copepod-2012_cmass_SON_gridded.mat','zoo_g')
cope_SON = zoo_g;
clear zoo_g

%% check orientation
figure
pcolor(cope_all)
shading flat
title('all')

figure
pcolor(cope_DJF)
shading flat
title('winter')

figure
pcolor(cope_MAM)
shading flat
title('spring')

figure
pcolor(cope_JJA)
shading flat
title('summer')

figure 
pcolor(cope_SON) 
shading flat
title('fall')

figure
pcolor(deptho)
shading flat
title('depth')

%% 
close all

%% Units
% From mg m-3 to mg m-2 (integrate top 200m)
dep200 = min(200,deptho);
dep200(isnan(deptho(:))) = nan;

cope_all = cope_all .* dep200;
cope_DJF = cope_DJF .* dep200;
cope_MAM = cope_MAM .* dep200;
cope_JJA = cope_JJA .* dep200;
cope_SON = cope_SON .* dep200;

%% CMIP6 models
load('cmip6_hist_space_means_50yr_seasons_zmeso200_glmm100_same_orientation.mat');
load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat')

%% check orientations
figure(1)
pcolor(gzmo_MAM); shading flat;
title('gfdl')

figure(6)
pcolor(lat_g); shading flat;
title('lat g')

%% 
close all

%% Vectorize, put in ABC order
% convert to same units
%models mol C/m2 -> g C/m2
%copepod mg C/m2 -> g C/m2
mod_all(:,1) = (czmo_all(:)) * 12.01;
mod_all(:,2) = (mzmo_all(:)) * 12.01;
mod_all(:,3) = (nzmo_all(:)) * 12.01;
mod_all(:,4) = (gzmo_all(:)) * 12.01;
mod_all(:,5) = (izmo_all(:)) * 12.01;
mod_all(:,6) = (uzmo_all(:)) * 12.01;
mod_all(:,7) = (cope_all(:)) * 1e-3;

mod_DJF(:,1) = (czmo_DJF(:)) * 12.01;
mod_DJF(:,2) = (mzmo_DJF(:)) * 12.01;
mod_DJF(:,3) = (nzmo_DJF(:)) * 12.01;
mod_DJF(:,4) = (gzmo_DJF(:)) * 12.01;
mod_DJF(:,5) = (izmo_DJF(:)) * 12.01;
mod_DJF(:,6) = (uzmo_DJF(:)) * 12.01;
mod_DJF(:,7) = (cope_DJF(:)) * 1e-3;

mod_JJA(:,1) = (czmo_JJA(:)) * 12.01;
mod_JJA(:,2) = (mzmo_JJA(:)) * 12.01;
mod_JJA(:,3) = (nzmo_JJA(:)) * 12.01;
mod_JJA(:,4) = (gzmo_JJA(:)) * 12.01;
mod_JJA(:,5) = (izmo_JJA(:)) * 12.01;
mod_JJA(:,6) = (uzmo_JJA(:)) * 12.01;
mod_JJA(:,7) = (cope_JJA(:)) * 1e-3;

mod_MAM(:,1) = (czmo_MAM(:)) * 12.01;
mod_MAM(:,2) = (mzmo_MAM(:)) * 12.01;
mod_MAM(:,3) = (nzmo_MAM(:)) * 12.01;
mod_MAM(:,4) = (gzmo_MAM(:)) * 12.01;
mod_MAM(:,5) = (izmo_MAM(:)) * 12.01;
mod_MAM(:,6) = (uzmo_MAM(:)) * 12.01;
mod_MAM(:,7) = (cope_MAM(:)) * 1e-3;

mod_SON(:,1) = (czmo_SON(:)) * 12.01;
mod_SON(:,2) = (mzmo_SON(:)) * 12.01;
mod_SON(:,3) = (nzmo_SON(:)) * 12.01;
mod_SON(:,4) = (gzmo_SON(:)) * 12.01;
mod_SON(:,5) = (izmo_SON(:)) * 12.01;
mod_SON(:,6) = (uzmo_SON(:)) * 12.01;
mod_SON(:,7) = (cope_SON(:)) * 1e-3;


%% all clim
lat = lat_g(:);
lon = lon_g(:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:9) = mod_all;

% nn = ~isnan(comb(:,8));
% comb = comb(nn,:);

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','copepod'});
writetable(obsmod,[sfile 'skill_hist_model_copepod_fullyr_clim_200.csv'])

%% Winter
lat = lat_g(:);
lon = lon_g(:);
dcomb(:,1) = lat;
dcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
dcomb(nh,3:9) = mod_DJF(nh,:);
dcomb(sh,3:9) = mod_JJA(sh,:);

nn = ~isnan(dcomb(:,8));
dcomb = dcomb(nn,:);

omDJF = array2table(dcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','copepod'});
writetable(omDJF,[sfile 'skill_hist_model_copepod_DJF_clim_200.csv'])

%% Summer clim
lat = lat_g(:);
lon = lon_g(:);

jcomb(:,1) = lat;
jcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
jcomb(nh,3:9) = mod_JJA(nh,:);
jcomb(sh,3:9) = mod_DJF(sh,:);

nn = ~isnan(jcomb(:,8));
jcomb = jcomb(nn,:);

omJJA = array2table(jcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','copepod'});
writetable(omJJA,[sfile 'skill_hist_model_copepod_JJA_clim_200.csv'])

%% Spring clim
lat = lat_g(:);
lon = lon_g(:);

mcomb(:,1) = lat;
mcomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
mcomb(nh,3:9) = mod_MAM(nh,:);
mcomb(sh,3:9) = mod_SON(sh,:);

nn = ~isnan(mcomb(:,8));
mcomb = mcomb(nn,:);

omMAM = array2table(mcomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','copepod'});
writetable(omMAM,[sfile 'skill_hist_model_copepod_MAM_clim_200.csv'])

%% Fall clim
lat = lat_g(:);
lon = lon_g(:);

scomb(:,1) = lat;
scomb(:,2) = lon;

nh = (lat>=0);
sh = (lat<0);
scomb(nh,3:9) = mod_SON(nh,:);
scomb(sh,3:9) = mod_MAM(sh,:);

nn = ~isnan(scomb(:,8));
scomb = scomb(nn,:);

omSON = array2table(scomb,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','copepod'});
writetable(omSON,[sfile 'skill_hist_model_copepod_SON_clim_200.csv'])

%%
save('skill_hist_model_copepod_climatols_seasons.mat','comb','dcomb','jcomb','mcomb','scomb')
