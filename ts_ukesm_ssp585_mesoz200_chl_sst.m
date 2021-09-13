% Save ts 
% UKESM mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'])
load([hpath 'ukesm_isimip_hist_sst_monthly_1950_2014.mat'])
load([hpath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'])

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
load([spath 'ukesm_isimip_ssp585_zmeso_200_monthly_2015_2100.mat'])
load([spath 'ukesm_isimip_ssp585_sst_monthly_2015_2100.mat'])
load([spath 'ukesm_isimip_ssp585_surf_chl_monthly_2015_2100.mat'])

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
HUzm = nanmean(reshape(Hzm_mo,ni*nj,nyh));
HUsst = nanmean(reshape(Hsst_mo,ni*nj,nyh));
HUchl = nanmean(reshape(Hchl_mo,ni*nj,nyh));

FUzm = nanmean(reshape(Fzm_mo,ni*nj,nyf));
FUsst = nanmean(reshape(Fsst_mo,ni*nj,nyf));
FUchl = nanmean(reshape(Fchl_mo,ni*nj,nyf));

%% Plots
figure(1)
subplot(2,2,1)
plot(Hyr(6:12:end),HUzm,'b'); hold on;
plot(Fyr(6:12:end),FUzm,'b');
title('mesoz')

subplot(2,2,3)
plot(Hyr(6:12:end),HUchl,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FUchl,'color',[0 0.75 0.5]);
title('schl')

subplot(2,2,4)
plot(Hyr(6:12:end),HUsst,'r'); hold on;
plot(Fyr(6:12:end),FUsst,'r');
title('sst')
print('-dpng',[ppath 'ukesm_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%%
save([fpath 'ukesm_hist_ssp585_tsmeans_zmeso200_schl_sst.mat'],'Hyr','Fyr',...
    'HUzm','HUsst','HUchl','FUzm','FUsst','FUchl');
