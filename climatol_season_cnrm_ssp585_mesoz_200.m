% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% Use regridded output
% Only use last 50 yrs (2051-2100)

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';

load([fpath 'cnrm_ssp585_zmeso200_monthly_onedeg_2015_2100.mat'])

%%
zmeso200(zmeso200>=1e19) = NaN;
ztest = zmeso200(:,:,10);
figure
pcolor(ztest)

%% Choose last 50 yrs to be similar time scale as historic
yid = find(time>2051);

zm50 = zmeso200(:,:,yid);

%% By hand
nmo = length(yid);
nyr = length(yid)/12;

zm_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zm50(:,:,mo)),3);
end

%% Seasonal clim
zmo_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmo_MAM = mean(zm_mo(:,:,3:5),3);
zmo_JJA = mean(zm_mo(:,:,6:8),3);
zmo_SON = mean(zm_mo(:,:,9:11),3);
zmo_all = mean(zm_mo,3);

%%
save([fpath 'cnrm_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],'yr',...
    'long_name','standard_name','units','lat','lon','yid',...
    'zmo_DJF','zmo_MAM','zmo_JJA','zmo_SON','zmo_all','zm_mo');
