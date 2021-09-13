% Save ts 
% UKESM npp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'ukesm_isimip_hist_npp_monthly_1951_2014.mat'])

%%
npp(npp>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hnpp = (npp);

%% annual means by grid cell
[ni,nj,ht] = size(Hnpp);

nyh = ht/12;

st = 1:12:ht;
en = 12:12:ht;

Hnpp_mo = nan*ones(ni,nj,nyh);
for m = 1:nyh
    mo = st(m):en(m);
    Hnpp_mo(:,:,m) = nanmean(Hnpp(:,:,mo),3);
end

%% SSP 
load([spath 'ukesm_isimip_ssp585_npp_monthly_2015_2100.mat'])

%% 
npp(npp>=1e19) = NaN;

yr = ((time)/12)+1601;
Fyr = yr;
Fnpp = (npp);

%% annual means by grid cell
[fi,fj,ft] = size(Fnpp);

nyf = ft/12;

st = 1:12:ft;
en = 12:12:ft;

Fnpp_mo = nan*ones(fi,fj,nyf);
for m = 1:nyf
    mo = st(m):en(m);
    Fnpp_mo(:,:,m) = nanmean(Fnpp(:,:,mo),3);
end

%% time means
HUnpp = nanmean(reshape(Hnpp_mo,ni*nj,nyh));

FUnpp = nanmean(reshape(Fnpp_mo,fi*fj,nyf));

%% Plots
figure(1)
plot(Hyr(6:12:end),HUnpp,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FUnpp,'color',[0 0.75 0.5]);
title('npp')
print('-dpng',[ppath 'ukesm_hist_ssp585_tsmeans_npp.png'])

%%
save([fpath 'ukesm_hist_ssp585_tsmeans_npp.mat'],'Hyr','Fyr',...
    'HUnpp','FUnpp');
