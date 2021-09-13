% Calc seasonal climatology
% DON'T Use ClimateDataToolbox

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

load([fpath 'cnrm_hist_zmeso_vint_monthly_1950_2014.mat'])

%% Monthly clim
yr = [yr1;yr2];
%runs = find(yr>1950 & yr<=2015);
runs = 1:length(yr);
[zmc,tc] = climatology(zmeso_vint,runs,'monthly');

%% By hand
nmo = length(runs);
nyr = length(runs)/12;

zm_mo = nan*ones(362,294,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zmeso_vint(:,:,mo)),3);
end

%% Seasonal clim
zmc_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmc_MAM = mean(zm_mo(:,:,3:5),3);
zmc_JJA = mean(zm_mo(:,:,6:8),3);
zmc_SON = mean(zm_mo(:,:,9:11),3);
zmc_all = mean(zm_mo,3);

%%
save([fpath 'cnrm_hist_zmeso_vint_climatol_1950_2014.mat'],'yr',...
    'long_name','standard_name','units_vint','lat','lon','runs',...
    'zmc_DJF','zmc_MAM','zmc_JJA','zmc_SON','zmc_all');
