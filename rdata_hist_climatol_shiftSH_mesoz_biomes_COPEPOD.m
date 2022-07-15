% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

%% depth info
Cdir = '/Volumes/MIP/Fish-MIP/CMIP6/';
load([Cdir 'GFDL/gridspec_gfdl_cmip6.mat'],'deptho','LAT','LON','lmask');

%% COPEPOD ---------------------------------------
load('copepod-2012_cmass_all_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

%% climatologies
cope_DJF = nanmean(zoo_g(:,:,[1 2 12]),3);
cope_MAM = nanmean(zoo_g(:,:,3:5),3);
cope_JJA = nanmean(zoo_g(:,:,6:8),3);
cope_SON = nanmean(zoo_g(:,:,9:11),3);
cope_all = nanmean(zoo_g,3);

% %% GLOBCOLOR ---------------------------------------
% % CHL
% opath ='/Volumes/MIP/Obs_Data/Chl/';
% load([opath 'globcolour_soblend_mc_x1.mat'])
%
% %% Reshape
% [lat_c,lon_c] = meshgrid(lat,lon);
%
% % climatologies
% schl_all = nanmean(chl,3);
%
% %% BIOME

%% MODISAQUA
% CHL
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load('glm100_obs_chl.mat');
load('glm100_obs_grid.mat')

% Reshape
[ni,nj] = size(lon_g);
glmc_mo = reshape(glm100obschl,ni,nj,12);
lat_mc = reshape(Lat,ni,nj);
lon_mc = reshape(Lon,ni,nj);

% climatologies
mchl_all = nanmean(glmc_mo,3);

%% BIOME
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';
load([bpath 'data_biomes_MODISAqua_x1.mat']);
[lat_mb,lon_mb] = meshgrid(lat,lon);

mbiomes = biomes;

clear lat lon biomes

%% SEAWIFS ---------------------------------------
% CHL
cpath ='/Volumes/MIP/Obs_Data/Chl/';
load([cpath 'SeaWiFS_Mission_Climatology.mat'])

%% Reshape
[lat_sc,lon_sc] = meshgrid(lat,lon);

% climatologies
schl_all = chl_all;

clear lat lon chl_all

%% BIOME
bpath2='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';

load([bpath2 'SeaWiFS_based_biomes_x1.mat']);
[lat_sb,lon_sb] = meshgrid(lat,lon);

clear lat lon

%% check orientations
figure(1)
pcolor(zoo_g); shading flat;
title('zoo')

figure(2)
pcolor(schl_all); shading flat;
title('Schl')

figure(3)
pcolor(mchl_all); shading flat;
title('Mchl')

figure(4)
pcolor(sbiomes); shading flat;
title('sb')

figure(5)
pcolor(mbiomes); shading flat;
title('mb')

figure(6)
pcolor(lat_g); shading flat;
title('lat g')

figure(7)
pcolor(lat_mc); shading flat;
title('lat mc')

figure(8)
pcolor(lat_mb); shading flat;
title('lat mb')

figure(9)
pcolor(lat_sc); shading flat;
title('lat sc')

figure(10)
pcolor(lat_sb); shading flat;
title('lat sb')

figure(11)
pcolor(deptho); shading flat;
title('depth')

figure(12)
pcolor(LAT); shading flat;
title('LAT')

%% Flip
close all

mbiomes = fliplr(mbiomes);
lat_g = fliplr(lat_g);
lat_mc = fliplr(lat_mc);
lat_mb = fliplr(lat_mb);

%% Units
% From mg m-3 to mg m-2 (integrate top 200m)
dep200 = min(200,deptho);
dep200(isnan(deptho(:))) = nan;
zoo_200 = zoo_g .* dep200;

%% save same orientation
save('space_means_copepod_zmeso200_chls_same_orientation.mat',...
    'zoo_200','schl_all','mchl_all','sbiomes','mbiomes','lat_sc','lon_sc');

%% Vectorize,
comb(:,1) = lat_sc(:);
comb(:,2) = lon_sc(:);
comb(:,3) = zoo_200(:) .* 1e-3; % From mg to g
comb(:,4) = schl_all(:);
comb(:,5) = mchl_all(:);
comb(:,6) = sbiomes(:);
comb(:,7) = mbiomes(:);

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','zmeso200','SEAWIFSchl','MODISchl','SEAWIFSbiomes','MODISbiomes'});
writetable(obsmod,'skill_hist_COPEPOD_all_clim_200_chls.csv')

%%
cope_comb = comb;
save('skill_hist_obs_copepod_data.mat','cope_comb')
