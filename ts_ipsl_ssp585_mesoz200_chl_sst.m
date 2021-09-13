% Save ts 
% IPSL mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'])
load([hpath 'ipsl_hist_sst_monthly_1950_2014.mat'])
load([hpath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'])

%%
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hzm = double(zmeso_200);
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
clear zmeso_200 sst schl
load([spath 'ipsl_ssp585_zmeso_200_monthly_2015_2100.mat'])
load([spath 'ipsl_ssp585_sst_monthly_2015_2100.mat'])
load([spath 'ipsl_ssp585_surf_chl_monthly_2015_2100.mat'])

%% 
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

Fyr = yr;
Fzm = double(zmeso_200);
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
HIzm = nanmean(reshape(Hzm_mo,ni*nj,nyh));
HIsst = nanmean(reshape(Hsst_mo,ni*nj,nyh));
HIchl = nanmean(reshape(Hchl_mo,ni*nj,nyh));

FIzm = nanmean(reshape(Fzm_mo,ni*nj,nyf));
FIsst = nanmean(reshape(Fsst_mo,ni*nj,nyf));
FIchl = nanmean(reshape(Fchl_mo,ni*nj,nyf));

%% Plots 
figure(1)
subplot(2,2,1)
plot(Hyr(6:12:end),HIzm,'b'); hold on;
plot(Fyr(6:12:end),FIzm,'b');
title('mesoz')

subplot(2,2,3)
plot(Hyr(6:12:end),HIchl,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FIchl,'color',[0 0.75 0.5]);
title('schl')

subplot(2,2,4)
plot(Hyr(6:12:end),HIsst,'r'); hold on;
plot(Fyr(6:12:end),FIsst,'r');
title('sst')
print('-dpng',[ppath 'ipsl_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%%
save([fpath 'ipsl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat'],'Hyr','Fyr',...
    'HIzm','HIsst','HIchl','FIzm','FIsst','FIchl');
