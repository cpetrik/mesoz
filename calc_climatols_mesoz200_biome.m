% CMIP6 output 
% 200m integrations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

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

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo','units');

cz = zm_mo;
cz(cz(:)<0) = 0;
cunits = units;

clear zm_mo units

%% CMCC
load([mpath 'cmcc_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo','units');

mz = zm_mo;
mz(mz(:)<0) = 0;
munits = units;

clear zm_mo units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo','units');

nz = zm_mo;
nz(nz(:)<0) = 0;
nunits = units;

clear zm_mo units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo');

uz = zm_mo;
uz(uz(:)<0) = 0;

clear zm_mo units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo');

iz = zm_mo;
iz(iz(:)<0) = 0;

clear zm_mo units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zm_mo','lat','lon');

gz = zm_mo;
gz(gz(:)<0) = 0;

[lat_g,lon_g] = meshgrid(lat,lon);

clear zm_mo 

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_mesoz.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj,nt] = size(gz);
glmz_mo = reshape(glmobsmesoz,ni,nj,12); 
lat_c = reshape(Lat,ni,nj);
lon_c = reshape(Lon,ni,nj);

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
cz = cz * 12.01 * 1e3;
mz = mz * 12.01 * 1e3;
nz = nz * 12.01 * 1e3;
gz = gz * 12.01 * 1e3;
iz = iz * 12.01 * 1e3;
uz = uz * 12.01 * 1e3;
%obsglm mg C/m2 

%% Check all same grids
figure %Ant R
pcolor(squeeze(cz(:,:,1)))
shading flat
title('CAN')

figure %Rus R
pcolor(hcbiomes)
shading flat
title('CAN')

figure %Ant R
pcolor(squeeze(mz(:,:,1)))
shading flat
title('Cm')

figure %Rus R
pcolor(hmbiomes)
shading flat
title('Cm')

figure %Ant R
pcolor(squeeze(nz(:,:,1)))
shading flat
title('CNRM')

figure %Rus R
pcolor(hnbiomes)
shading flat
title('CNRM')

figure %Rus R
pcolor(squeeze(gz(:,:,1)))
shading flat
title('GFDL')

figure %Rus R
pcolor(hgbiomes)
shading flat
title('GFDL')

figure %Ant R
pcolor(squeeze(iz(:,:,1)))
shading flat
title('IPSL')

figure %Rus R
pcolor(hibiomes)
shading flat
title('IPSL')

figure %Ant R
pcolor(squeeze(uz(:,:,1)))
shading flat
title('UK')

figure %Rus R
pcolor(hubiomes)
shading flat
title('UK')

figure %Rus R
pcolor(squeeze(glmz_mo(:,:,6)))
shading flat
title('obs')

figure %Ant R
pcolor(biomes)
shading flat
title('obs')


%% Flip as needed
cz = fliplr(cz);
mz = fliplr(mz);
nz = fliplr(nz);
iz = fliplr(iz);
uz = fliplr(uz);
biomes = fliplr(biomes);

%% grids
close all
figure %Rus R
pcolor(lat_b)
shading flat
title('lat b')

figure %Rus R
pcolor(lat_c)
shading flat
title('lat c')

figure %Ant R
pcolor(lat_m)
shading flat
title('lat m')

figure %Rus R
pcolor(lat_g)
shading flat
title('lat g')

%% Save fixed
save('Hist_mesoz_clim_biomes_samegrid.mat','cz','mz','nz','gz','iz','uz','glmz_mo',...
    'hcbiomes','hmbiomes','hnbiomes','hgbiomes','hibiomes','hubiomes','biomes',...
    'lat_g','lon_g')

%% Reshape
cz = reshape(cz,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
nz = reshape(nz,ni*nj,nt);
gz = reshape(gz,ni*nj,nt);
iz = reshape(iz,ni*nj,nt);
uz = reshape(uz,ni*nj,nt);
oz = reshape(glmz_mo,ni*nj,nt);

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

%% Take means by biome
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

%%
Cz = nanmean(Scz);
Mz = nanmean(Smz);
Nz = nanmean(Snz);
Gz = nanmean(Sgz);
Iz = nanmean(Siz);
Uz = nanmean(Suz);
Oz = nanmean(Soz);

Czlc = nanmean(Scz(cid,:));
Mzlc = nanmean(Smz(mid,:));
Nzlc = nanmean(Snz(nid,:));
Gzlc = nanmean(Sgz(gid,:));
Izlc = nanmean(Siz(iid,:));
Uzlc = nanmean(Suz(uid,:));
Ozlc = nanmean(Soz(oid,:));

Czss = nanmean(Scz(css,:));
Mzss = nanmean(Smz(mss,:));
Nzss = nanmean(Snz(nss,:));
Gzss = nanmean(Sgz(gss,:));
Izss = nanmean(Siz(iss,:));
Uzss = nanmean(Suz(uss,:));
Ozss = nanmean(Soz(oss,:));

Czps = nanmean(Scz(cps,:));
Mzps = nanmean(Smz(mps,:));
Nzps = nanmean(Snz(nps,:));
Gzps = nanmean(Sgz(gps,:));
Izps = nanmean(Siz(ips,:));
Uzps = nanmean(Suz(ups,:));
Ozps = nanmean(Soz(ops,:));

%%
save('Hist_ts_mesoz_clim_biomes.mat','Cz','Mz','Nz','Gz','Iz','Uz','Oz',...
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
ylim([0 1250])










