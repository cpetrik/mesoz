% Correlations of ESM and obs chl
% obs chl in mg/m3 and model in kg/m3?

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_surf_chl_monthly_onedeg_1951_2014.mat'],'schl',...
    'units','runs','yr');

cz = nanmean(schl,3) * 1e3;
cz(cz(:)<0) = 0;
cunits = units;
ctime = yr(runs);

clear schl units runs yr

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_surf_chl_monthly_onedeg_1951_2014.mat'],...
    'schl','units','runs','yr');

nz = nanmean(schl,3);
nz(nz(:)<0) = 0;
nunits = units;
ntime = yr(runs);

clear schl units runs yr

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],...
    'schl','units','runs','yr');

uz = nanmean(schl,3) * 1e3;
uz(uz(:)<0) = 0;
uunits = units;
utime = yr(runs);

clear schl units runs yr

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat']);

iz = nanmean(schl,3);
iz(iz(:)<0) = 0;
iunits = units;
itime = yr(runs);
[lat_g,lon_g] = meshgrid(lat,lon);

clear schl units runs yr lat lon
 
%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],...
    'schl','units','runs','yr');

gz = nanmean(schl,3) * 1e3;
gz(gz(:)<0) = 0;
gunits = units;
gtime = yr(runs);

gz = fliplr(gz);
pcolor(gz);

clear schl units runs yr 

%% Chl
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';
opath='/Volumes/MIP/Obs_Data/Chl/';
load([opath 'globcolour_soblend_mc_x1.mat'])

chl_ocx = chl*1e-3;
[lat_c,lon_c] = meshgrid(lat,lon);

mchl = nanmean(chl_ocx,3);

figure
pcolor(lon_c,lat_c,mchl);
figure
pcolor(mchl);

%% Scatter plots
%add corr line

x=-6:0.25:-1;

figure(4)
subplot(2,3,1)
scatter(log10(mchl(:)),log10(cz(:))); hold on
plot(x,x,'--k')
axis([-6 -1 -6 -1])

subplot(2,3,2)
scatter(log10(mchl(:)),log10(nz(:)));hold on
plot(x,x,'--k')
axis([-6 -1 -6 -1])

subplot(2,3,3)
scatter(log10(mchl(:)),log10(gz(:)));hold on
plot(x,x,'--k')
axis([-6 -1 -6 -1])

subplot(2,3,4)
scatter(log10(mchl(:)),log10(iz(:)));hold on
plot(x,x,'--k')
axis([-6 -1 -6 -1])

subplot(2,3,5)
scatter(log10(mchl(:)),log10(uz(:)));hold on
plot(x,x,'--k')
axis([-6 -1 -6 -1])


%% Vectorize, put in ABC order
model(:,1) = log10(cz(:));
model(:,2) = log10(nz(:));
model(:,3) = log10(gz(:));
model(:,4) = log10(iz(:));
model(:,5) = log10(uz(:));
obs = log10(mchl(:));

%%
nn = ~isnan(obs);
model = model(nn,:);
obs = obs(nn);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3) = obs;
comb(:,4:8) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','obs','CAN','CNRM','GFDL','IPSL','UK'});

save('chl_obs_mod_mean_1950_2014_all_lats.mat','model','obs','obsmod');
writetable(obsmod,'chl_obs_mod_mean_1950_2014_all_lats.csv')





