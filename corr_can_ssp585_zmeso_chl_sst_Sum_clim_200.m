% Correlations of zmeso biomass with surf chl and sst
% CAN-OE

clear all
close all

%% CAN
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';
load([npath 'can_ssp585_zmeso_200_climatol_2051_2100_gridded.mat']);

%%
zJ(zJ(:)<0) = 0;
zunits = units;

%%
load([npath 'can_ssp585_chl_sst_climatol_2051_2100_gridded.mat'],...
    'cJ','sJ','units');

cunits = units;

%% Vectorize
model(:,1) = (zJ(:));
model(:,2) = (cJ(:));
model(:,3) = (sJ(:));

%%
nn = ~isnan(zM);
model = model(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:5) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','zooJJA','chlJJA','sstJJA'});

save('can_ssp585_mod_zoo_chl_sst_Sum_clim_200.mat','model','comb','obsmod');
writetable(obsmod,'can_ssp585_mod_zoo_chl_sst_Sum_clim_200.csv')

%% Scatter plots
figure(1)
scatter(log10(model(:,2)),log10(model(:,1)),'k'); hold on










