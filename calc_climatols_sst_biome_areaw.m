% CMIP6 output 
% SST
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';

load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

%% Biomes
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

load([fpath 'CanESM_historical_biomes.mat'],'hcbiomes');
load([mpath 'CMCC_historical_biomes.mat'],'hmbiomes');
load([fpath 'CNRM_historical_biomes.mat'],'hnbiomes');
load([fpath 'GFDL_historical_biomes.mat'],'hgbiomes');
load([fpath 'IPSL_historical_biomes.mat'],'hibiomes');
load([fpath 'UKESM_historical_biomes.mat'],'hubiomes','lat','lon');

[lat_b,lon_b] = meshgrid(lat,lon);
clear lat lon

load([bpath 'data_biomes_MODISAqua_x1.mat']);
[lat_m,lon_m] = meshgrid(lat,lon);

%% Load fixed
load('Hist_sst_clim_biomes_samegrid.mat')

%% Reshape
[ni,nj,nm] = size(gt);
ct = reshape(ct,ni*nj,nm);
mt = reshape(mt,ni*nj,nm);
nt = reshape(nt,ni*nj,nm);
gt = reshape(gt,ni*nj,nm);
it = reshape(it,ni*nj,nm);
ut = reshape(ut,ni*nj,nm);
ot = reshape(ot_mo,ni*nj,nm);

area = repmat(cell_area,1,1,nm);
areav = reshape(area,ni*nj,nm);

%% Shift SH by 6mo
nh = (lat_g(:)>=0);
sh = (lat_g(:)<0);

Sct(nh,:)    = ct(nh,:);
Sct(sh,1:6)  = ct(sh,7:12);
Sct(sh,7:12) = ct(sh,1:6);

Smt(nh,:)    = mt(nh,:);
Smt(sh,1:6)  = mt(sh,7:12);
Smt(sh,7:12) = mt(sh,1:6);

Snt(nh,:)    = nt(nh,:);
Snt(sh,1:6)  = nt(sh,7:12);
Snt(sh,7:12) = nt(sh,1:6);

Sgt(nh,:)    = gt(nh,:);
Sgt(sh,1:6)  = gt(sh,7:12);
Sgt(sh,7:12) = gt(sh,1:6);

Sit(nh,:)    = it(nh,:);
Sit(sh,1:6)  = it(sh,7:12);
Sit(sh,7:12) = it(sh,1:6);

Sut(nh,:)    = ut(nh,:);
Sut(sh,1:6)  = ut(sh,7:12);
Sut(sh,7:12) = ut(sh,1:6);

Sot(nh,:)    = ot(nh,:);
Sot(sh,1:6)  = ot(sh,7:12);
Sot(sh,7:12) = ot(sh,1:6);

%% ID each biome
% Exclude LC > 45N/S
pol = find(lat_g <= 45 & lat_g >= -45);

clc = find(hcbiomes==1);
mlc = find(hmbiomes==1);
nlc = find(hnbiomes==1);
glc = find(hgbiomes==1);
ilc = find(hibiomes==1);
ulc = find(hubiomes==1);
olc = find(biomes==1);

cid = intersect(clc,pol);
mid = intersect(mlc,pol);
nid = intersect(nlc,pol);
gid = intersect(glc,pol);
iid = intersect(ilc,pol);
uid = intersect(ulc,pol);
oid = intersect(olc,pol);

css = (hcbiomes==2);
mss = (hmbiomes==2);
nss = (hnbiomes==2);
gss = (hgbiomes==2);
iss = (hibiomes==2);
uss = (hubiomes==2);
oss = (biomes==2);

cps = (hcbiomes==3);
mps = (hmbiomes==3);
nps = (hnbiomes==3);
gps = (hgbiomes==3);
ips = (hibiomes==3);
ups = (hubiomes==3);
ops = (biomes==3);

%% Take means by biome
% Ct = nanmean(Sct);
% Mt = nanmean(Smt);
% Nt = nanmean(Snt);
% Gt = nanmean(Sgt);
% It = nanmean(Sit);
% Ut = nanmean(Sut);
% Ot = nanmean(Sot);

Ct = nansum(Sct .* areav) ./ nansum(areav);
Mt = nansum(Smt .* areav) ./ nansum(areav);
Nt = nansum(Snt .* areav) ./ nansum(areav);
Gt = nansum(Sgt .* areav) ./ nansum(areav);
It = nansum(Sit .* areav) ./ nansum(areav);
Ut = nansum(Sut .* areav) ./ nansum(areav);
Ot = nansum(Sot .* areav) ./ nansum(areav);

Ctlc = nansum(Sct(cid,:) .* areav(cid,:)) ./ nansum(areav(cid,:));
Mtlc = nansum(Smt(mid,:) .* areav(mid,:)) ./ nansum(areav(mid,:));
Ntlc = nansum(Snt(nid,:) .* areav(nid,:)) ./ nansum(areav(nid,:));
Gtlc = nansum(Sgt(gid,:) .* areav(gid,:)) ./ nansum(areav(gid,:));
Itlc = nansum(Sit(iid,:) .* areav(iid,:)) ./ nansum(areav(iid,:));
Utlc = nansum(Sut(uid,:) .* areav(uid,:)) ./ nansum(areav(uid,:));
Otlc = nansum(Sot(oid,:) .* areav(oid,:)) ./ nansum(areav(oid,:));

Ctss = nansum(Sct(css,:) .* areav(css,:)) ./ nansum(areav(css,:));
Mtss = nansum(Smt(mss,:) .* areav(mss,:)) ./ nansum(areav(mss,:));
Ntss = nansum(Snt(nss,:) .* areav(nss,:)) ./ nansum(areav(nss,:));
Gtss = nansum(Sgt(gss,:) .* areav(gss,:)) ./ nansum(areav(gss,:));
Itss = nansum(Sit(iss,:) .* areav(iss,:)) ./ nansum(areav(iss,:));
Utss = nansum(Sut(uss,:) .* areav(uss,:)) ./ nansum(areav(uss,:));
Otss = nansum(Sot(oss,:) .* areav(oss,:)) ./ nansum(areav(oss,:));

Ctps = nansum(Sct(cps,:) .* areav(cps,:)) ./ nansum(areav(cps,:));
Mtps = nansum(Smt(mps,:) .* areav(mps,:)) ./ nansum(areav(mps,:));
Ntps = nansum(Snt(nps,:) .* areav(nps,:)) ./ nansum(areav(nps,:));
Gtps = nansum(Sgt(gps,:) .* areav(gps,:)) ./ nansum(areav(gps,:));
Itps = nansum(Sit(ips,:) .* areav(ips,:)) ./ nansum(areav(ips,:));
Utps = nansum(Sut(ups,:) .* areav(ups,:)) ./ nansum(areav(ups,:));
Otps = nansum(Sot(ops,:) .* areav(ops,:)) ./ nansum(areav(ops,:));

%%
save('Hist_ts_sst_clim_biomes_areaw.mat','Ct','Mt','Nt','Gt','It','Ut','Ot',...
    'Ctlc','Mtlc','Ntlc','Gtlc','Itlc','Utlc','Otlc',...
    'Ctss','Mtss','Ntss','Gtss','Itss','Utss','Otss',...
    'Ctps','Mtps','Ntps','Gtps','Itps','Utps','Otps')

%% Plots 
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey

cb=[34/255 136/255 51/255;...   %green
    %238/255 119/255 51/255;...  %orange
    %0/255 153/255 136/255;...   %teal
    153/255 153/255 51/255;...   %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];                     %black

set(groot,'defaultAxesColorOrder',cb);

%%
mo = 1:12;

close all

figure(1)
subplot(4,3,1)
plot(mo,Ct); hold on;
plot(mo,Mt); hold on;
plot(mo,Nt); hold on;
plot(mo,Gt); hold on;
plot(mo,It); hold on;
plot(mo,Ut); hold on;
plot(mo,Ot); hold on;
title('SST (^oC)')
ylabel('Global')
xlim([1 12])

subplot(4,3,4)
plot(mo,Ctlc); hold on;
plot(mo,Mtlc); hold on;
plot(mo,Ntlc); hold on;
plot(mo,Gtlc); hold on;
plot(mo,Itlc); hold on;
plot(mo,Utlc); hold on;
plot(mo,Otlc); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,7)
plot(mo,Ctss); hold on;
plot(mo,Mtss); hold on;
plot(mo,Ntss); hold on;
plot(mo,Gtss); hold on;
plot(mo,Itss); hold on;
plot(mo,Utss); hold on;
plot(mo,Otss); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,10)
plot(mo,Ctps); hold on;
plot(mo,Mtps); hold on;
plot(mo,Ntps); hold on;
plot(mo,Gtps); hold on;
plot(mo,Itps); hold on;
plot(mo,Utps); hold on;
plot(mo,Otps); hold on;
ylabel('HCPS')
xlim([1 12])









