% Create vectors of zmeso biomass with surf chl and sst
% Labeled with biomes for use in R
% UKESM

clear all
close all

%% UKESM biome_masks
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';
load([fpath 'UKESM_historical_biomes.mat'],'biomes');

%% UKESM zmeso
npath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([npath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat']);
% Climatol of last 50 yrs
zmo_all(zmo_all(:)<0) = 0;
zunits = units;

%% UKESM chl and sst
load([npath 'ukesm_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all','units');

cunits = units;

%% check orientation
figure(1)
pcolor(zmo_all); shading flat;

figure(2)
pcolor(chl_all); shading flat;

figure(3)
pcolor(sst_all); shading flat;

figure(4)
pcolor(biomes); shading flat;

[lat_g,lon_g] = meshgrid(lat,lon);

figure(5)
pcolor(lat_g); shading flat;

figure(6)
pcolor(lon_g); shading flat;

%% Vectorize
biomes = fliplr(biomes);
sst_all = fliplr(sst_all);

model(:,1) = (biomes(:));
model(:,2) = (zmo_all(:));
model(:,3) = (chl_all(:));
model(:,4) = (sst_all(:));

%%
nn = ~isnan(zmo_all);
model = model(nn,:);

%%
lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:6) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','biomes','zmeso','chl','sst'});

save('ukesm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.mat','model','comb','obsmod');
writetable(obsmod,'ukesm_hist_zmeso200_chl_sst_clim_1965_2014_biomes.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,3)),log10(model(:,2)),10,model(:,1),'filled'); hold on
