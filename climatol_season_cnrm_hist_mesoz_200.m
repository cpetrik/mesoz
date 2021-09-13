% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% Use regridded output
% Only use last 50 yrs (1965-2014)

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';

load([fpath 'cnrm_hist_zmeso200_monthly_onedeg_1951_2014.mat'])

%% Monthly clim
%runs = find(time>1965);
runs = 181:length(yr);
zmeso50 = zmeso200(:,:,runs);

%% By hand
nmo = length(runs);
nyr = length(runs)/12;

zm_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zmeso50(:,:,mo)),3);
end

%% Seasonal clim
zmo_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmo_MAM = mean(zm_mo(:,:,3:5),3);
zmo_JJA = mean(zm_mo(:,:,6:8),3);
zmo_SON = mean(zm_mo(:,:,9:11),3);
zmo_all = mean(zm_mo,3);

%%
save([fpath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],'yr',...
    'long_name','standard_name','units','lat','lon','runs',...
    'zmo_DJF','zmo_MAM','zmo_JJA','zmo_SON','zmo_all','zm_mo');
