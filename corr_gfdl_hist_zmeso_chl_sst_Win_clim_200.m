% Correlations of zmeso biomass with surf chl and sst
% GFDL

clear all
close all

%% GFDL
npath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([npath 'gfdl_hist_zmeso_200_climatol_1950_2014_gridded.mat']);

%
zD(zD(:)<0) = 0;
zunits = units;

%
load([npath 'gfdl_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'cD','sD','units');

cunits = units;

%% Vectorize
model(:,1) = (zD(:));
model(:,2) = (cD(:));
model(:,3) = (sD(:));

%%
nn = ~isnan(zD);
model = model(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:5) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','zooDJF','chlDJF','sstDJF'});

save('gfdl_hist_mod_zoo_chl_sst_Win_clim_200.mat','model','comb','obsmod');
writetable(obsmod,'gfdl_hist_mod_zoo_chl_sst_Win_clim_200.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,2)),log10(model(:,1)),'k'); hold on










