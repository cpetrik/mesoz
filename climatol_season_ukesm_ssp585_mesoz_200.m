% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% UKESM mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

load([fpath 'ukesm_isimip_ssp585_zmeso_200_monthly_2015_2100.mat'])

%%
zmeso_200(zmeso_200>=1e19) = NaN;

ztest = zmeso_200(:,:,10);

figure
pcolor(ztest)

%% Choose last 50 yrs to be similar time scale as historic
yid = find(yr>2051);

zm50 = zmeso_200(:,:,yid);

%% By hand
nmo = length(yid);
nyr = length(yid)/12;

zm_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zm50(:,:,mo)),3);
end

%% By hand shows seasonal progression between hemispheres
zmo_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmo_MAM = mean(zm_mo(:,:,3:5),3);
zmo_JJA = mean(zm_mo(:,:,6:8),3);
zmo_SON = mean(zm_mo(:,:,9:11),3);
zmo_all = mean(zm_mo,3);

%%
save([fpath 'ukesm_ssp585_zmeso200_onedeg_climatol_2051_2100.mat'],'yr',...
    'long_name','standard_name','units_vint','lat','lon','yid',...
    'zmo_DJF','zmo_MAM','zmo_JJA','zmo_SON','zmo_all','zm_mo');
