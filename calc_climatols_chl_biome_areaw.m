% CMIP6 output 
% surf chl
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';

load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

%% Biomes
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';
mpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

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
load('Hist_chl_clim_biomes_samegrid.mat')

%% Reshape
[ni,nj,nt] = size(gc);
cc = reshape(cc,ni*nj,nt);
mc = reshape(mc,ni*nj,nt);
nc = reshape(nc,ni*nj,nt);
gc = reshape(gc,ni*nj,nt);
ic = reshape(ic,ni*nj,nt);
uc = reshape(uc,ni*nj,nt);
oc = reshape(oc_mo,ni*nj,nt);

area = repmat(cell_area,1,1,nt);
areav = reshape(area,ni*nj,nt);

%% Shift SH by 6mo
nh = (lat_g(:)>=0);
sh = (lat_g(:)<0);

Scc(nh,:)    = cc(nh,:);
Scc(sh,1:6)  = cc(sh,7:12);
Scc(sh,7:12) = cc(sh,1:6);

Smc(nh,:)    = mc(nh,:);
Smc(sh,1:6)  = mc(sh,7:12);
Smc(sh,7:12) = mc(sh,1:6);

Snc(nh,:)    = nc(nh,:);
Snc(sh,1:6)  = nc(sh,7:12);
Snc(sh,7:12) = nc(sh,1:6);

Sgc(nh,:)    = gc(nh,:);
Sgc(sh,1:6)  = gc(sh,7:12);
Sgc(sh,7:12) = gc(sh,1:6);

Sic(nh,:)    = ic(nh,:);
Sic(sh,1:6)  = ic(sh,7:12);
Sic(sh,7:12) = ic(sh,1:6);

Suc(nh,:)    = uc(nh,:);
Suc(sh,1:6)  = uc(sh,7:12);
Suc(sh,7:12) = uc(sh,1:6);

Soc(nh,:)    = oc(nh,:);
Soc(sh,1:6)  = oc(sh,7:12);
Soc(sh,7:12) = oc(sh,1:6);

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
% Cc = nanmean(Scc);
% Mc = nanmean(Smc);
% Nc = nanmean(Snc);
% Gc = nanmean(Sgc);
% Ic = nanmean(Sic);
% Uc = nanmean(Suc);
% Oc = nanmean(Soc);

Cc = nansum(Scc .* areav) ./ nansum(areav);
Mc = nansum(Smc .* areav) ./ nansum(areav);
Nc = nansum(Snc .* areav) ./ nansum(areav);
Gc = nansum(Sgc .* areav) ./ nansum(areav);
Ic = nansum(Sic .* areav) ./ nansum(areav);
Uc = nansum(Suc .* areav) ./ nansum(areav);
Oc = nansum(Soc .* areav) ./ nansum(areav);

Cclc = nansum(Scc(cid,:) .* areav(cid,:)) ./ nansum(areav(cid,:));
Mclc = nansum(Smc(mid,:) .* areav(mid,:)) ./ nansum(areav(mid,:));
Nclc = nansum(Snc(nid,:) .* areav(nid,:)) ./ nansum(areav(nid,:));
Gclc = nansum(Sgc(gid,:) .* areav(gid,:)) ./ nansum(areav(gid,:));
Iclc = nansum(Sic(iid,:) .* areav(iid,:)) ./ nansum(areav(iid,:));
Uclc = nansum(Suc(uid,:) .* areav(uid,:)) ./ nansum(areav(uid,:));
Oclc = nansum(Soc(oid,:) .* areav(oid,:)) ./ nansum(areav(oid,:));

Ccss = nansum(Scc(css,:) .* areav(css,:)) ./ nansum(areav(css,:));
Mcss = nansum(Smc(mss,:) .* areav(mss,:)) ./ nansum(areav(mss,:));
Ncss = nansum(Snc(nss,:) .* areav(nss,:)) ./ nansum(areav(nss,:));
Gcss = nansum(Sgc(gss,:) .* areav(gss,:)) ./ nansum(areav(gss,:));
Icss = nansum(Sic(iss,:) .* areav(iss,:)) ./ nansum(areav(iss,:));
Ucss = nansum(Suc(uss,:) .* areav(uss,:)) ./ nansum(areav(uss,:));
Ocss = nansum(Soc(oss,:) .* areav(oss,:)) ./ nansum(areav(oss,:));

Ccps = nansum(Scc(cps,:) .* areav(cps,:)) ./ nansum(areav(cps,:));
Mcps = nansum(Smc(mps,:) .* areav(mps,:)) ./ nansum(areav(mps,:));
Ncps = nansum(Snc(nps,:) .* areav(nps,:)) ./ nansum(areav(nps,:));
Gcps = nansum(Sgc(gps,:) .* areav(gps,:)) ./ nansum(areav(gps,:));
Icps = nansum(Sic(ips,:) .* areav(ips,:)) ./ nansum(areav(ips,:));
Ucps = nansum(Suc(ups,:) .* areav(ups,:)) ./ nansum(areav(ups,:));
Ocps = nansum(Soc(ops,:) .* areav(ops,:)) ./ nansum(areav(ops,:));

%%
save('Hist_ts_chl_clim_biomes_areaw.mat','Cc','Mc','Nc','Gc','Ic','Uc','Oc',...
    'Cclc','Mclc','Nclc','Gclc','Iclc','Uclc','Oclc',...
    'Ccss','Mcss','Ncss','Gcss','Icss','Ucss','Ocss',...
    'Ccps','Mcps','Ncps','Gcps','Icps','Ucps','Ocps')

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
plot(mo,Cc); hold on;
plot(mo,Mc); hold on;
plot(mo,Nc); hold on;
plot(mo,Gc); hold on;
plot(mo,Ic); hold on;
plot(mo,Uc); hold on;
plot(mo,Oc); hold on;
title('chl (mg m^-^2)')
ylabel('Global')
xlim([1 12])

subplot(4,3,4)
plot(mo,Cclc); hold on;
plot(mo,Mclc); hold on;
plot(mo,Nclc); hold on;
plot(mo,Gclc); hold on;
plot(mo,Iclc); hold on;
plot(mo,Uclc); hold on;
plot(mo,Oclc); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,7)
plot(mo,Ccss); hold on;
plot(mo,Mcss); hold on;
plot(mo,Ncss); hold on;
plot(mo,Gcss); hold on;
plot(mo,Icss); hold on;
plot(mo,Ucss); hold on;
plot(mo,Ocss); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,10)
plot(mo,Ccps); hold on;
plot(mo,Mcps); hold on;
plot(mo,Ncps); hold on;
plot(mo,Gcps); hold on;
plot(mo,Icps); hold on;
plot(mo,Ucps); hold on;
plot(mo,Ocps); hold on;
ylabel('HCPS')
xlim([1 12])









