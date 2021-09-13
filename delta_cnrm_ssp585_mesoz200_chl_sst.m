% Save time chunks to calc change in future 
% CNRM mesoz regridded by J Luo

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'cnrm_hist_zmeso200_monthly_onedeg_1951_2014.mat'])
load([hpath 'cnrm_hist_sst_monthly_onedeg_1951_2014.mat'])
load([hpath 'cnrm_hist_surf_chl_monthly_onedeg_1951_2014.mat'])

%%
zmeso200(zmeso200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hzm = double(zmeso200);
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
load([spath 'cnrm_ssp585_zmeso200_monthly_onedeg_2015_2100.mat'])
load([spath 'cnrm_ssp585_sst_monthly_onedeg_2015_2100.mat'])
load([spath 'cnrm_ssp585_surf_chl_monthly_onedeg_2015_2100.mat'])

%% 
zmeso200(zmeso200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

Fyr = yr;
Fzm = double(zmeso200);
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

%% space means
hyr = Hyr(6:12:end);
fyr = Fyr(6:12:end);
% last 50 yrs: 1965-2014 and 2051-2100
h50 = find(hyr>=1965 & hyr<2015);
f50 = find(fyr>=2051 & fyr<2101);

%%
HNzm50 = nanmean(Hzm_mo(:,:,h50),3);
HNsst50 = nanmean(Hsst_mo(:,:,h50),3);
HNchl50 = nanmean(Hchl_mo(:,:,h50),3);

FNzm50 = nanmean(Fzm_mo(:,:,f50),3);
FNsst50 = nanmean(Fsst_mo(:,:,f50),3);
FNchl50 = nanmean(Fchl_mo(:,:,f50),3);

%% first 10 yrs: 1965-1974 and last 10 yrs: 2091-2100
h10 = find(hyr>=1965 & hyr<1975);
f10 = find(fyr>=2091 & fyr<2101);

HNzm10 = nanmean(Hzm_mo(:,:,h10),3);
HNsst10 = nanmean(Hsst_mo(:,:,h10),3);
HNchl10 = nanmean(Hchl_mo(:,:,h10),3);

FNzm10 = nanmean(Fzm_mo(:,:,f10),3);
FNsst10 = nanmean(Fsst_mo(:,:,f10),3);
FNchl10 = nanmean(Fchl_mo(:,:,f10),3);

%% 90s: 1991-2000 and last 10 yrs: 2091-2100
h90 = find(hyr>=1991 & hyr<2001);
f90 = find(fyr>=2091 & fyr<2101);

HNzm90 = nanmean(Hzm_mo(:,:,h90),3);
HNsst90 = nanmean(Hsst_mo(:,:,h90),3);
HNchl90 = nanmean(Hchl_mo(:,:,h90),3);

FNzm90 = nanmean(Fzm_mo(:,:,f90),3);
FNsst90 = nanmean(Fsst_mo(:,:,f90),3);
FNchl90 = nanmean(Fchl_mo(:,:,f90),3);

%%
save([fpath 'cnrm_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'hyr','fyr',...
    'h50','f50','h10','f10','HNzm50','HNsst50','HNchl50',...
    'FNzm50','FNsst50','FNchl50','HNzm10','HNsst10','HNchl10',...
    'FNzm10','FNsst10','FNchl10','HNzm90','HNsst90','HNchl90',...
    'FNzm90','FNsst90','FNchl90');
