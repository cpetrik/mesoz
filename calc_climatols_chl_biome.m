% CMIP6 output 
% surf chl

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

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

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo');
cc = chl_mo;
cc(cc(:)<0) = 0;
clear chl_mo 

%% CMCC
load([mpath 'cmcc_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo');
mc = chl_mo;
mc(mc(:)<0) = 0;
clear chl_mo 

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo');
nc = chl_mo;
nc(nc(:)<0) = 0;
clear chl_mo 

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo');
uc = chl_mo;
uc(uc(:)<0) = 0;
clear chl_mo 

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo');
ic = chl_mo;
ic(ic(:)<0) = 0;
clear chl_mo 

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_mo','lat','lon');
gc = chl_mo;
gc(gc(:)<0) = 0;
[lat_g,lon_g] = meshgrid(lat,lon);
clear chl_mo 

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_chl.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj,nt] = size(gc);
oc_mo = reshape(glmobschl,ni,nj,12); 
lat_c = reshape(Lat,ni,nj);
lon_c = reshape(Lon,ni,nj);

%% Convert units 
%chl in kg m-3', put in mg m-3 (except CNRM & IPSL already that)
cc = cc * 1e6;
mc = mc * 1e6;
gc = gc * 1e6;
uc = uc * 1e6;
nc = nc * 1e3;
ic = ic * 1e3;

%% Check all same grids
figure %Ant R
pcolor(squeeze(cc(:,:,1)))
shading flat
title('CAN')

figure %Rus R
pcolor(hcbiomes)
shading flat
title('CAN')

figure %Ant R
pcolor(squeeze(mc(:,:,1)))
shading flat
title('CMCC')

figure %Rus R
pcolor(hmbiomes)
shading flat
title('CMCC')

figure %Ant R
pcolor(squeeze(nc(:,:,1)))
shading flat
title('CNRM')

figure %Rus R
pcolor(hnbiomes)
shading flat
title('CNRM')

figure %Rus R
pcolor(squeeze(gc(:,:,1)))
shading flat
title('GFDL')

figure %Rus R
pcolor(hgbiomes)
shading flat
title('GFDL')

figure %Ant R
pcolor(squeeze(ic(:,:,1)))
shading flat
title('IPSL')

figure %Rus R
pcolor(hibiomes)
shading flat
title('IPSL')

figure %Ant R
pcolor(squeeze(uc(:,:,1)))
shading flat
title('UK')

figure %Rus R
pcolor(hubiomes)
shading flat
title('UK')

figure %Rus R
pcolor(squeeze(oc_mo(:,:,6)))
shading flat
title('obs')

figure %Ant R
pcolor(biomes)
shading flat
title('obs')


%% Flip as needed
cc = fliplr(cc);
mc = fliplr(mc);
nc = fliplr(nc);
ic = fliplr(ic);
uc = fliplr(uc);
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
save('Hist_chl_clim_biomes_samegrid.mat','cc','mc','nc','gc','ic','uc','oc_mo',...
    'hcbiomes','hmbiomes','hnbiomes','hgbiomes','hibiomes','hubiomes','biomes',...
    'lat_g','lon_g')


%% Reshape
cc = reshape(cc,ni*nj,nt);
mc = reshape(mc,ni*nj,nt);
nc = reshape(nc,ni*nj,nt);
gc = reshape(gc,ni*nj,nt);
ic = reshape(ic,ni*nj,nt);
uc = reshape(uc,ni*nj,nt);
oc = reshape(oc_mo,ni*nj,nt);

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
Cc = nanmean(Scc);
Mc = nanmean(Smc);
Nc = nanmean(Snc);
Gc = nanmean(Sgc);
Ic = nanmean(Sic);
Uc = nanmean(Suc);
Oc = nanmean(Soc);

Cclc = nanmean(Scc(cid,:));
Mclc = nanmean(Smc(mid,:));
Nclc = nanmean(Snc(nid,:));
Gclc = nanmean(Sgc(gid,:));
Iclc = nanmean(Sic(iid,:));
Uclc = nanmean(Suc(uid,:));
Oclc = nanmean(Soc(oid,:));

Ccss = nanmean(Scc(css,:));
Mcss = nanmean(Smc(mss,:));
Ncss = nanmean(Snc(nss,:));
Gcss = nanmean(Sgc(gss,:));
Icss = nanmean(Sic(iss,:));
Ucss = nanmean(Suc(uss,:));
Ocss = nanmean(Soc(oss,:));

Ccps = nanmean(Scc(cps,:));
Mcps = nanmean(Smc(mps,:));
Ncps = nanmean(Snc(nps,:));
Gcps = nanmean(Sgc(gps,:));
Icps = nanmean(Sic(ips,:));
Ucps = nanmean(Suc(ups,:));
Ocps = nanmean(Soc(ops,:));

%%
save('Hist_ts_chl_clim_biomes.mat','Cc','Mc','Nc','Gc','Ic','Uc','Oc',...
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









