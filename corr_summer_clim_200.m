% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zJ','units');

cz = zJ;
cz(cz(:)<0) = 0;
cunits = units;

clear zJ units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zJ','units');

nz = zJ;
nz(nz(:)<0) = 0;
nunits = units;

clear zJ units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zJ','units');

uz = zJ;
uz(uz(:)<0) = 0;
uunits = units;

clear zJ units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zJ','units');

iz = zJ;
iz(iz(:)<0) = 0;
iunits = units;

clear zJ units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zJ','units');

gz = zJ;
gz(gz(:)<0) = 0;
gunits = units;

clear zJ units

%% COPEPOD
load('copepod-2012_cmass_JJA_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

% From m-3 to m-2 (integrate top 200m)
zoo_200 = zoo_g*200;

%% Chl
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';

load([fpath 'seawifs_chl_ocx_growingseason_mean.mat'])

[lat_c,lon_c] = meshgrid(lat,lon);

chl_ocx = fliplr(chl_ocx);
lat_c = fliplr(lat_c);
lon_c = fliplr(lon_c);

chl = chl_ocx(:,11:end);

%% Vectorize, put in ABC order
model(:,1) = log10(cz(:));
model(:,2) = log10(nz(:));
model(:,3) = log10(gz(:));
model(:,4) = log10(iz(:));
model(:,5) = log10(uz(:));
obs = log10(zoo_200(:));
chl2 = log10(chl(:));

%%
nn = ~isnan(obs);
model = model(nn,:);
obs = obs(nn);
chl2 = chl2(nn);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3) = obs;
comb(:,4:8) = model;
comb(:,9) = chl2;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','obs','CAN','CNRM','GFDL','IPSL','UK','chl'});

save('obs_mod_chl_summer_clim_200.mat','model','obs','obsmod');
writetable(obsmod,'obs_mod_chl_summer_clim_200.csv')

%% Scatter plots
%add corr line

x=-10:0.1:5;

figure(4)
subplot(2,3,1)
scatter(obs,model(:,1)); hold on
plot(x,x,'--k')
axis([0 5 -10 4.5])

subplot(2,3,2)
scatter(obs,model(:,2));hold on
plot(x,x,'--k')
axis([0 5 1.5 4.5])

subplot(2,3,3)
scatter(obs,model(:,3));hold on
plot(x,x,'--k')
axis([0.5 5 1.5 4])

subplot(2,3,4)
scatter(obs,model(:,4));hold on
plot(x,x,'--k')
axis([0 5 1 3.5])

subplot(2,3,5)
scatter(obs,model(:,5));hold on
plot(x,x,'--k')
axis([0 5 0 4])

subplot(2,3,6)
scatter(obs,chl2);hold on
axis([0 5 -2 1])

