% Save ts 
% IPSL npp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'ipsl_hist_npp_monthly_1951_2014.mat'])

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
load([spath 'ipsl_ssp585_npp_monthly_2015_2100.mat'])

%% 
npp(npp>=1e19) = NaN;

yr = ((time)/12)+1601;
Fyr = yr;
Fnpp = (npp);

%% annual means by grid cell
[ni,nj,ft] = size(Fnpp);

nyf = ft/12;

st = 1:12:ft;
en = 12:12:ft;

Fnpp_mo = nan*ones(ni,nj,nyf);
for m = 1:nyf
    mo = st(m):en(m);
    Fnpp_mo(:,:,m) = nanmean(Fnpp(:,:,mo),3);
end

%% time means
HInpp = nanmean(reshape(Hnpp_mo,ni*nj,nyh));

FInpp = nanmean(reshape(Fnpp_mo,ni*nj,nyf));

%% Plots
figure(1)
plot(Hyr(6:12:end),HInpp,'color',[0 0.75 0.5]); hold on;
plot(Fyr(6:12:end),FInpp,'color',[0 0.75 0.5]);
title('npp')
print('-dpng',[ppath 'ipsl_hist_ssp585_tsmeans_npp.png'])

%%
save([fpath 'ipsl_hist_ssp585_tsmeans_npp.mat'],'Hyr','Fyr',...
    'HInpp','FInpp');
