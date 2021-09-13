% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% Use regridded output
% Only use last 50 yrs (2051-2100)

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/ssp585/';

load([fpath 'cmcc_ssp585_sst_monthly_onedeg_2015_2100.mat'])
load([fpath 'cmcc_ssp585_surf_chl_monthly_onedeg_2015_2100.mat'])

%%
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

ttest = sst(:,:,20);
ctest = schl(:,:,20);

figure
pcolor(ttest)
figure
pcolor(ctest)

%% Choose last 50 yrs to be similar time scale as historic
yid = find(time>2051);

sst5 = sst(:,:,yid);
chl5 = schl(:,:,yid);

%% By hand
nmo = length(yid);
nyr = length(yid)/12;

sst_mo = nan*ones(360,180,12);
chl_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    sst_mo(:,:,m) = nanmean(double(sst5(:,:,mo)),3);
    chl_mo(:,:,m) = nanmean(double(chl5(:,:,mo)),3);
end

%% Seasonal clim
sst_DJF = mean(sst_mo(:,:,[1 2 12]),3);
sst_MAM = mean(sst_mo(:,:,3:5),3);
sst_JJA = mean(sst_mo(:,:,6:8),3);
sst_SON = mean(sst_mo(:,:,9:11),3);
sst_all = mean(sst_mo,3);

chl_DJF = mean(chl_mo(:,:,[1 2 12]),3);
chl_MAM = mean(chl_mo(:,:,3:5),3);
chl_JJA = mean(chl_mo(:,:,6:8),3);
chl_SON = mean(chl_mo(:,:,9:11),3);
chl_all = mean(chl_mo,3);

%%
save([fpath 'cmcc_ssp585_chl_sst_onedeg_climatol_2051_2100.mat'],'time',...
    'long_name','standard_name','units','lat','lon','yid',...
    'sst_DJF','sst_MAM','sst_JJA','sst_SON','sst_all','sst_mo',...
    'chl_DJF','chl_MAM','chl_JJA','chl_SON','chl_all','chl_mo');
