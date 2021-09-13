% Calc seasonal climatology
% DON'T Use ClimateDataToolbox
% Use regridded output
% Only use last 50 yrs (1965-2014)

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%%
load([fpath 'cmcc_hist_sst_monthly_onedeg_1965_2014.mat'])
load([fpath 'cmcc_hist_surf_chl_monthly_onedeg_1965_2014.mat'])

%% By hand
[ni,nj,nmo] = size(sst);
nyr = nmo/12;

sst_mo = nan*ones(ni,nj,12);
chl_mo = nan*ones(ni,nj,12);
for m = 1:12
    mo = m:12:nyr;
    sst_mo(:,:,m) = nanmean((sst(:,:,mo)),3);
    chl_mo(:,:,m) = nanmean((schl(:,:,mo)),3);
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
save([fpath 'cmcc_hist_chl_sst_onedeg_climatol_1965_2014.mat'],'time',...
    'long_name','standard_name','units','lat','lon',...
    'sst_DJF','sst_MAM','sst_JJA','sst_SON','sst_all','sst_mo',...
    'chl_DJF','chl_MAM','chl_JJA','chl_SON','chl_all','chl_mo');

