% Save ts 
% GFDL mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_sst_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'])

%%
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hzm = double(zmeso_200);
Hsst = double(sst);
Hschl = double(schl);

%% annual means by grid cell
[ni,nj,ht] = size(Hzm);

nyh = ht/12;

st = 1:12:ht;
en = 12:12:ht;

Hzm_mo = nan*ones(ni,nj,nyh);
Hsst_mo = nan*ones(ni,nj,nyh);
Hchl_mo = nan*ones(ni,nj,nyh);
for m = 1:nyh
    mo = st(m):en(m);
    Hzm_mo(:,:,m) = nanmean(Hzm(:,:,mo),3);
    Hsst_mo(:,:,m) = nanmean(Hsst(:,:,mo),3);
    Hchl_mo(:,:,m) = nanmean(Hschl(:,:,mo),3);
end

%% SSP 
load([spath 'gfdl_ssp585_zmeso_200_monthly_2015_2100.mat'])
load([spath 'gfdl_ssp585_sst_monthly_2015_2100.mat'])
load([spath 'gfdl_ssp585_surf_chl_monthly_2015_2100.mat'])

%% 
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

Fyr = yr;
Fzm = double(zmeso_200);
Fsst = double(sst);
Fschl = double(schl);

%% annual means by grid cell
[ni,nj,ft] = size(Fzm);

nyf = ft/12;

st = 1:12:ft;
en = 12:12:ft;

Fzm_mo = nan*ones(ni,nj,nyf);
Fsst_mo = nan*ones(ni,nj,nyf);
Fchl_mo = nan*ones(ni,nj,nyf);
for m = 1:nyf
    mo = st(m):en(m);
    Fzm_mo(:,:,m) = nanmean(Fzm(:,:,mo),3);
    Fsst_mo(:,:,m) = nanmean(Fsst(:,:,mo),3);
    Fchl_mo(:,:,m) = nanmean(Fschl(:,:,mo),3);
end

%% time means
HGzm = nanmean(reshape(Hzm_mo,ni*nj,nyh));
HGsst = nanmean(reshape(Hsst_mo,ni*nj,nyh));
HGchl = nanmean(reshape(Hchl_mo,ni*nj,nyh));

FGzm = nanmean(reshape(Fzm_mo,ni*nj,nyf));
FGsst = nanmean(reshape(Fsst_mo,ni*nj,nyf));
FGchl = nanmean(reshape(Fchl_mo,ni*nj,nyf));

%% Plots
figure(1)
subplot(2,2,1)
plot(Hyr(6:12:end),HGzm,'b'); hold on;
plot(Fyr(6:12:end),FGzm,'b');
title('mesoz')

subplot(2,2,3)
plot(Hyr(6:12:end),HGchl,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FGchl,'color',[0 0.75 0.5]);
title('schl')

subplot(2,2,4)
plot(Hyr(6:12:end),HGsst,'r'); hold on;
plot(Fyr(6:12:end),FGsst,'r');
title('sst')
print('-dpng',[ppath 'gfdl_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%%
save([fpath 'gfdl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat'],'Hyr','Fyr',...
    'HGzm','HGsst','HGchl','FGzm','FGsst','FGchl');
