% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

%% Stromberg obs model
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/Stromberg_x1_all/';
load([fpath 'StrombergQTR-m00_int200_mgCm2.mat']);
[lat_m,lon_m] = meshgrid(lat,lon);

szmo_all = zmo_all;

clear lat lon

%% Biomes
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';

load([bpath 'SeaWiFS_based_biomes_x1.mat']);
[lat_b,lon_b] = meshgrid(lat,lon);

%% check orientations
figure(1)
pcolor(szmo_all); shading flat;
title('Obs')

figure(2)
pcolor(sbiomes); shading flat;
title('obs')

figure(3)
pcolor(lat_m); shading flat;
title('Zlat')

figure(4)
pcolor(lat_b); shading flat;
title('Blat')

%% Flip
close all

szmo_all = fliplr(szmo_all);
lat_m = fliplr(lat_m);

%% save same orientation
save('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat',...
    'szmo_all','-append');

%% Vectorize, put in ABC order
% convert to same units
%models mol C/m2 -> g C/m2
%obsglm mg C/m2 -> g C/m2
mod_all = (szmo_all(:)) * 1e-3;

%% ADD BIOMES
bvec = sbiomes(:);

%% all clim
lat = lat_b(:);
lon = lon_b(:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3) = mod_all;
comb(:,4) = bvec;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','Stromberg','obsbiome'});
writetable(obsmod,'skill_hist_Stromberg_all_clim_200.csv')

%%
strom_comb = comb;
save('skill_hist_model_obsglm100_climatols.mat','strom_comb','-append')

