% Create vectors of zmeso biomass with surf chl and sst
% Labeled with biomes for use in R

% Chl and SST are those used in GLM to estimate mesoz
% satellite climatologies from MODIS-AQUA
% NOT sure if chl matches up with biomes

clear all
close all

%% biome_masks
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/';
load([fpath 'data_biomes_x1.mat']);
[lat_b,lon_b] = meshgrid(lat,lon);

clear lat lon

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_chl.mat']);
load([opath 'glm_obs_sst.mat']);
load([opath 'glm_obs_mesoz.mat'])
load([opath 'glm_obs_grid.mat'])

%% Reshape
[ni,nj] = size(biomes);
chl = reshape(glmobschl,ni,nj,12);
sst = reshape(glmobssst,ni,nj,12);
zmeso = reshape(glmobsmesoz,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);
bathy = reshape(Bathy,ni,nj);

%% take annual means 
mchl = nanmean(chl,3);
msst = nanmean(sst,3);
mzoo = nanmean(zmeso,3);

%% check orientation
figure(1)
pcolor(mzoo); shading flat;

figure(2)
pcolor(mchl); shading flat;

figure(3)
pcolor(msst); shading flat;

figure(4)
pcolor(biomes); shading flat;

%%
figure(5)
pcolor(lat_g); shading flat;
title('zoo lat')

figure(7)
pcolor(lat_b); shading flat;
title('biome lat')

%%
figure(6)
pcolor(lon_g); shading flat;
title('zoo lon')

figure(8)
pcolor(lon_b); shading flat;
title('biome lon')


%% Flip
biomes = fliplr(biomes);
lat_b = fliplr(lat_b);
lon_b = fliplr(lon_b);

%% Vectorize
model(:,1) = (biomes(:));
model(:,2) = (mzoo(:));
model(:,3) = (mchl(:));
model(:,4) = (msst(:));

%%
%nn = ~isnan(mzoo);
nn = ~isnan(biomes);
model = model(nn,:);

lat = lat_b(:);
lon = lon_b(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:6) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','biomes','zmeso','chl','sst'});

save('glm_obs_zmeso_chl_sst_ann_mean_biomes.mat','model','comb','obsmod');
writetable(obsmod,'glm_obs_zmeso_chl_sst_ann_mean_biomes.csv')

%% Scatter plots
figure(10)
scatter(log10(model(:,3)),log10(model(:,2)),10,model(:,1),'filled'); hold on
