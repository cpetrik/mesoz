% Save ts 
% CAN mesoz regridded by J Luo

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'can_hist_zmeso200_monthly_onedeg_1951_2014.mat'])
load([hpath 'can_hist_sst_monthly_onedeg_1951_2014.mat'])
load([hpath 'can_hist_surf_chl_monthly_onedeg_1951_2014.mat'])

%%
zmeso200(zmeso200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hzm = double(zmeso200);
Hsst = double(sst);
Hschl = double(schl);

%%
figure
pcolor(squeeze(Hzm(:,:,10)))
shading flat

figure
pcolor(squeeze(Hsst(:,:,10)))
shading flat

figure
pcolor(squeeze(Hschl(:,:,10)))
shading flat

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
load([spath 'can_ssp585_zmeso200_monthly_onedeg_2015_2100.mat'])
load([spath 'can_ssp585_sst_monthly_onedeg_2015_2100.mat'])
load([spath 'can_ssp585_surf_chl_monthly_onedeg_2015-2100.mat'])

%% 
zmeso200(zmeso200==0) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

Fyr = yr;
Fzm = double(zmeso200);
Fsst = double(sst);
Fschl = double(schl);

%%
figure
pcolor(squeeze(Fzm(:,:,15)))
shading flat

figure
pcolor(squeeze(Fsst(:,:,15)))
shading flat

figure
pcolor(squeeze(Fschl(:,:,15)))
shading flat

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
HCzm = nanmean(reshape(Hzm_mo,ni*nj,nyh));
HCsst = nanmean(reshape(Hsst_mo,ni*nj,nyh));
HCchl = nanmean(reshape(Hchl_mo,ni*nj,nyh));

FCzm = nanmean(reshape(Fzm_mo,ni*nj,nyf));
FCsst = nanmean(reshape(Fsst_mo,ni*nj,nyf));
FCchl = nanmean(reshape(Fchl_mo,ni*nj,nyf));

%% Plots 
figure(1)
subplot(2,2,1)
plot(Hyr(6:12:end),HCzm,'b'); hold on;
plot(Fyr(6:12:end),FCzm,'b');
title('mesoz')

subplot(2,2,3)
plot(Hyr(6:12:end),HCchl,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FCchl,'color',[0 0.75 0.5]);
title('schl')

subplot(2,2,4)
plot(Hyr(6:12:end),HCsst,'r'); hold on;
plot(Fyr(6:12:end),FCsst,'r');
title('sst')
print('-dpng',[ppath 'can_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%%
save([fpath 'can_hist_ssp585_tsmeans_zmeso200_schl_sst.mat'],'Hyr','Fyr',...
    'HCzm','HCsst','HCchl','FCzm','FCsst','FCchl');
