% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% Use regridded output
% Only use last 50 yrs (1965-2014)

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

load([fpath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'])

%% Monthly clim
time = yr(runs);
yid = find(time>1965);
zmeso50 = zmeso_200(:,:,yid);

%% By hand
nmo = length(yid);
nyr = length(yid)/12;

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
save([fpath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],'yr','time',...
    'long_name','standard_name','units_vint','lat','lon','runs','yid',...
    'zmo_DJF','zmo_MAM','zmo_JJA','zmo_SON','zmo_all','zm_mo');
