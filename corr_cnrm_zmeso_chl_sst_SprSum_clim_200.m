% Correlations of zmeso biomass with surf chl and sst
% CNRM

clear all
close all

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso_200_climatol_1950_2014_gridded.mat']);

%%
zM(zM(:)<0) = 0;
zJ(zJ(:)<0) = 0;
zunits = units;

%%
load([npath 'cnrm_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'cM','cJ','sM','sJ','units');

cunits = units;

%% Vectorize
model(:,1) = (zM(:));
model(:,2) = (cM(:));
model(:,3) = (sM(:));
model(:,4) = (zJ(:));
model(:,5) = (cJ(:));
model(:,6) = (sJ(:));

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
    {'Lat','Lon','zooMAM','chlMAM','sstMAM','zooJJA','chlJJA','sstJJA'});

save('cnrm_mod_zoo_chl_sst_SprSum_clim_200.mat','model','comb','obsmod');
writetable(obsmod,'cnrm_mod_zoo_chl_sst_SprSum_clim_200.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,2)),log10(model(:,1)),'k'); hold on

