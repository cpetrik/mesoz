% Create vectors of zmeso biomass with surf chl and sst
% Labeled with biomes for use in R

% Chl is 1997-2018 SeaWiFS and MODIS product + 
% the Johnson et al. Southern Ocean modification
% from J Luo
% Also used to construct biomes

% SST is 1971-2000 from OISST downloaded 2/16/21

clear all
close all

%% biome_masks
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/';
load([fpath 'data_biomes_x1.mat']);
[lat_b,lon_b] = meshgrid(lat,lon);

clear lat lon

% in vector form
vbiome = reshape(biomes,360*180,1);
vlat = reshape(lat_b,360*180,1);
vlon = reshape(lon_b,360*180,1);

%% Chl
cpath = '/Volumes/MIP/Obs_Data/Chl/';
load([cpath 'globcolour_soblend_mc_x1.mat']);
[lat_c,lon_c] = meshgrid(lat,lon);

clear lat lon

%% SST
spath = '/Volumes/MIP/Obs_Data/OISST/';
load([spath 'sst_annual_mean_onedeg_1971_2000.mat'],'Lat1','Lon1',...
    'sstG');

%% COPEPOD zmeso
load('copepod-2012_cmass_all_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

% From m-3 to m-2 (integrate top 200m)
zoo_200 = zoo_g*200;

%% take annual mean of chl
mchl = nanmean(chl,3);

%% check orientation
figure(1)
pcolor(zoo_200); shading flat;

figure(2)
pcolor(mchl); shading flat;

figure(3)
pcolor(sstG); shading flat;

figure(4)
pcolor(biomes); shading flat;

%%
figure(5)
pcolor(lat_g); shading flat;
title('zoo lat')

figure(7)
pcolor(lat_b); shading flat;
title('biome lat')

figure(9)
pcolor(lat_c); shading flat;
title('chl lat')

figure(11)
pcolor(Lat1); shading flat;
title('sst lat')

%%
figure(6)
pcolor(lon_g); shading flat;
title('zoo lon')

figure(8)
pcolor(lon_b); shading flat;
title('biome lon')

figure(10)
pcolor(lon_c); shading flat;
title('chl lon')

figure(12)
pcolor(Lon1); shading flat;
title('sst lon')

%% Flip
zoo_200 = fliplr(zoo_200);
sstG = fliplr(sstG);
lat_g = fliplr(lat_g);
Lat1 = fliplr(Lat1);
lon_g = fliplr(lon_g);
Lon1 = fliplr(Lon1);

%% Vectorize
model(:,1) = (biomes(:));
model(:,2) = (zoo_200(:));
model(:,3) = (mchl(:));
model(:,4) = (sstG(:));

%%
nn = ~isnan(zoo_200);
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

save('obs_zmeso200_chl_sst_clim_biomes.mat','model','comb','obsmod');
writetable(obsmod,'obs_zmeso200_chl_sst_clim_biomes.csv')

%% Scatter plots
figure(10)
scatter(log10(model(:,3)),log10(model(:,2)),10,model(:,1),'filled'); hold on
