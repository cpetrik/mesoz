% Save time chunks to calc change in future 
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

%% space means
hyr = Hyr(6:12:end);
fyr = Fyr(6:12:end);
% last 50 yrs: 1965-2014 and 2051-2100
h50 = find(hyr>=1965 & hyr<2015);
f50 = find(fyr>=2051 & fyr<2101);

%%
HGzm50 = nanmean(Hzm_mo(:,:,h50),3);
HGsst50 = nanmean(Hsst_mo(:,:,h50),3);
HGchl50 = nanmean(Hchl_mo(:,:,h50),3);

FGzm50 = nanmean(Fzm_mo(:,:,f50),3);
FGsst50 = nanmean(Fsst_mo(:,:,f50),3);
FGchl50 = nanmean(Fchl_mo(:,:,f50),3);

%% first 10 yrs: 1965-1974 and last 10 yrs: 2091-2100
h10 = find(hyr>=1965 & hyr<1975);
f10 = find(fyr>=2091 & fyr<2101);

HGzm10 = nanmean(Hzm_mo(:,:,h10),3);
HGsst10 = nanmean(Hsst_mo(:,:,h10),3);
HGchl10 = nanmean(Hchl_mo(:,:,h10),3);

FGzm10 = nanmean(Fzm_mo(:,:,f10),3);
FGsst10 = nanmean(Fsst_mo(:,:,f10),3);
FGchl10 = nanmean(Fchl_mo(:,:,f10),3);

%% 90s: 1991-2000 and last 10 yrs: 2091-2100
h90 = find(hyr>=1991 & hyr<2001);
f90 = find(fyr>=2091 & fyr<2101);

HGzm90 = nanmean(Hzm_mo(:,:,h90),3);
HGsst90 = nanmean(Hsst_mo(:,:,h90),3);
HGchl90 = nanmean(Hchl_mo(:,:,h90),3);

FGzm90 = nanmean(Fzm_mo(:,:,f90),3);
FGsst90 = nanmean(Fsst_mo(:,:,f90),3);
FGchl90 = nanmean(Fchl_mo(:,:,f90),3);

%%
save([fpath 'gfdl_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'hyr','fyr',...
    'h50','f50','h10','f10','HGzm50','HGsst50','HGchl50',...
    'FGzm50','FGsst50','FGchl50','HGzm10','HGsst10','HGchl10',...
    'FGzm10','FGsst10','FGchl10','HGzm90','HGsst90','HGchl90',...
    'FGzm90','FGsst90','FGchl90');


