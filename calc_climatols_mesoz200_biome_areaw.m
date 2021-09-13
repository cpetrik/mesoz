% CMIP6 output 
% 200m integrations
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

%% load fixed
load('Hist_mesoz_clim_biomes_samegrid.mat')

%% Reshape
[ni,nj,nt] = size(gz);
cz = reshape(cz,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
nz = reshape(nz,ni*nj,nt);
gz = reshape(gz,ni*nj,nt);
iz = reshape(iz,ni*nj,nt);
uz = reshape(uz,ni*nj,nt);
oz = reshape(glmz_mo,ni*nj,nt);

area = repmat(cell_area,1,1,nt);
areav = reshape(area,ni*nj,nt);

%% Shift SH by 6mo
nh = (lat_g(:)>=0);
sh = (lat_g(:)<0);

Scz(nh,:)    = cz(nh,:);
Scz(sh,1:6)  = cz(sh,7:12);
Scz(sh,7:12) = cz(sh,1:6);

Smz(nh,:)    = mz(nh,:);
Smz(sh,1:6)  = mz(sh,7:12);
Smz(sh,7:12) = mz(sh,1:6);

Snz(nh,:)    = nz(nh,:);
Snz(sh,1:6)  = nz(sh,7:12);
Snz(sh,7:12) = nz(sh,1:6);

Sgz(nh,:)    = gz(nh,:);
Sgz(sh,1:6)  = gz(sh,7:12);
Sgz(sh,7:12) = gz(sh,1:6);

Siz(nh,:)    = iz(nh,:);
Siz(sh,1:6)  = iz(sh,7:12);
Siz(sh,7:12) = iz(sh,1:6);

Suz(nh,:)    = uz(nh,:);
Suz(sh,1:6)  = uz(sh,7:12);
Suz(sh,7:12) = uz(sh,1:6);

Soz(nh,:)    = oz(nh,:);
Soz(sh,1:6)  = oz(sh,7:12);
Soz(sh,7:12) = oz(sh,1:6);

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
Cz = nansum(Scz .* areav) ./ nansum(areav);
Mz = nansum(Smz .* areav) ./ nansum(areav);
Nz = nansum(Snz .* areav) ./ nansum(areav);
Gz = nansum(Sgz .* areav) ./ nansum(areav);
Iz = nansum(Siz .* areav) ./ nansum(areav);
Uz = nansum(Suz .* areav) ./ nansum(areav);
Oz = nansum(Soz .* areav) ./ nansum(areav);

Czlc = nansum(Scz(cid,:) .* areav(cid,:)) ./ nansum(areav(cid,:));
Mzlc = nansum(Smz(mid,:) .* areav(mid,:)) ./ nansum(areav(mid,:));
Nzlc = nansum(Snz(nid,:) .* areav(nid,:)) ./ nansum(areav(nid,:));
Gzlc = nansum(Sgz(gid,:) .* areav(gid,:)) ./ nansum(areav(gid,:));
Izlc = nansum(Siz(iid,:) .* areav(iid,:)) ./ nansum(areav(iid,:));
Uzlc = nansum(Suz(uid,:) .* areav(uid,:)) ./ nansum(areav(uid,:));
Ozlc = nansum(Soz(oid,:) .* areav(oid,:)) ./ nansum(areav(oid,:));

Czss = nansum(Scz(css,:) .* areav(css,:)) ./ nansum(areav(css,:));
Mzss = nansum(Smz(mss,:) .* areav(mss,:)) ./ nansum(areav(mss,:));
Nzss = nansum(Snz(nss,:) .* areav(nss,:)) ./ nansum(areav(nss,:));
Gzss = nansum(Sgz(gss,:) .* areav(gss,:)) ./ nansum(areav(gss,:));
Izss = nansum(Siz(iss,:) .* areav(iss,:)) ./ nansum(areav(iss,:));
Uzss = nansum(Suz(uss,:) .* areav(uss,:)) ./ nansum(areav(uss,:));
Ozss = nansum(Soz(oss,:) .* areav(oss,:)) ./ nansum(areav(oss,:));

Czps = nansum(Scz(cps,:) .* areav(cps,:)) ./ nansum(areav(cps,:));
Mzps = nansum(Smz(mps,:) .* areav(mps,:)) ./ nansum(areav(mps,:));
Nzps = nansum(Snz(nps,:) .* areav(nps,:)) ./ nansum(areav(nps,:));
Gzps = nansum(Sgz(gps,:) .* areav(gps,:)) ./ nansum(areav(gps,:));
Izps = nansum(Siz(ips,:) .* areav(ips,:)) ./ nansum(areav(ips,:));
Uzps = nansum(Suz(ups,:) .* areav(ups,:)) ./ nansum(areav(ups,:));
Ozps = nansum(Soz(ops,:) .* areav(ops,:)) ./ nansum(areav(ops,:));

%%
save('Hist_ts_mesoz_clim_biomes_areaw.mat','Cz','Mz','Nz','Gz','Iz','Uz','Oz',...
    'Czlc','Mzlc','Nzlc','Gzlc','Izlc','Uzlc','Ozlc',...
    'Czss','Mzss','Nzss','Gzss','Izss','Uzss','Ozss',...
    'Czps','Mzps','Nzps','Gzps','Izps','Uzps','Ozps')

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
plot(mo,Cz); hold on;
plot(mo,Mz); hold on;
plot(mo,Nz); hold on;
plot(mo,Gz); hold on;
plot(mo,Iz); hold on;
plot(mo,Uz); hold on;
title('mesozoo (mgC m^-^2)')
ylabel('Global')
xlim([1 12])
%ylim([0 1100])

subplot(4,3,4)
plot(mo,Czlc); hold on;
plot(mo,Mzlc); hold on;
plot(mo,Nzlc); hold on;
plot(mo,Gzlc); hold on;
plot(mo,Izlc); hold on;
plot(mo,Uzlc); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,7)
plot(mo,Czss); hold on;
plot(mo,Mzss); hold on;
plot(mo,Nzss); hold on;
plot(mo,Gzss); hold on;
plot(mo,Izss); hold on;
plot(mo,Uzss); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,10)
plot(mo,Czps); hold on;
plot(mo,Mzps); hold on;
plot(mo,Nzps); hold on;
plot(mo,Gzps); hold on;
plot(mo,Izps); hold on;
plot(mo,Uzps); hold on;
ylabel('HCPS')
xlim([1 12])
%ylim([0 1500])










