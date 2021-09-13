% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% UKESM mesoz_vint prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

load([fpath 'ukesm_isimip_hist_zmeso_vint_monthly_1950_2014.mat'])

%% Monthly clim
runs = find(yr>1950 & yr<=2015);
[zmc,tc] = climatology(zmeso_vint,runs,'monthly');

%% By hand
nmo = length(runs);
nyr = length(runs)/12;

zm_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zmeso_vint(:,:,mo)),3);
end

%% By hand shows seasonal progression between hemispheres
zmc_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmc_MAM = mean(zm_mo(:,:,3:5),3);
zmc_JJA = mean(zm_mo(:,:,6:8),3);
zmc_SON = mean(zm_mo(:,:,9:11),3);
zmc_all = mean(zm_mo,3);

%%
save([fpath 'ukesm_isimip_hist_zmeso_vint_climatol_1950_2014.mat'],'yr',...
    'long_name','standard_name','units','lat','lon','runs',...
    'zmc_DJF','zmc_MAM','zmc_JJA','zmc_SON','zmc_all');

