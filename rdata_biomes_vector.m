% Create vectors of biomes for use in R

clear all
close all

%% biome_masks
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/';
load([fpath 'data_biomes_x1.mat']);
[lat_b,lon_b] = meshgrid(lat,lon);

%% Flip
biomes = fliplr(biomes);
lat_b = fliplr(lat_b);
lon_b = fliplr(lon_b);

%% Vectorize
model(:,1) = (lat_b(:));
model(:,2) = (lon_b(:));
model(:,3) = (biomes(:));

%%
nn = find(~isnan(model(:,3)));

obsmod = array2table(model,'VariableNames',...
    {'Lat','Lon','biomes'});

%%
save('obs_biomes_vector.mat','model','obsmod');
writetable(obsmod,'obs_biomes_vector.csv')
