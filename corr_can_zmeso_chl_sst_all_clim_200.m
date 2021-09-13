% Correlations of zmeso biomass with surf chl and sst
% CAN-OE

clear all
close all

%% CAN
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([npath 'can_hist_zmeso_200_climatol_1950_2014_gridded.mat']);

%%
zall(zall(:)<0) = 0;
zunits = units;

%%
load([npath 'can_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'call','sall','units');

cunits = units;

%% Vectorize
model(:,1) = (zall(:));
model(:,2) = (call(:));
model(:,3) = (sall(:));

%%
nn = ~isnan(zall);
model = model(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:5) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','zoo','chl','sst'});

save('can_mod_zoo_chl_sst_all_clim_200.mat','model','comb','obsmod');
writetable(obsmod,'can_mod_zoo_chl_sst_all_clim_200.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,2)),log10(model(:,1)),'k'); hold on










