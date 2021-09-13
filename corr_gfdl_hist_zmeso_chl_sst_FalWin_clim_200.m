% Correlations of zmeso biomass with surf chl and sst
% GFDL

clear all
close all

%% GFDL
npath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([npath 'gfdl_hist_zmeso_200_climatol_1950_2014_gridded.mat']);

%
zS(zS(:)<0) = 0;
zD(zD(:)<0) = 0;
zunits = units;

%
load([npath 'gfdl_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'cS','cD','sS','sD','units');

cunits = units;

%% Vectorize
model(:,1) = (zS(:));
model(:,2) = (cS(:));
model(:,3) = (sS(:));
model(:,4) = (zD(:));
model(:,5) = (cD(:));
model(:,6) = (sD(:));

%%
nn = ~isnan(zM);
model = model(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:8) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','zooSON','chlSON','sstSON','zooDJF','chlDJF','sstDJF'});

save('gfdl_hist_mod_zoo_chl_sst_FalWin_clim_200.mat','model','comb','obsmod');
writetable(obsmod,'gfdl_hist_mod_zoo_chl_sst_FalWin_clim_200.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,2)),log10(model(:,1)),'k'); hold on










