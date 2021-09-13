% Maps of ESM and obs biomes
% Hist and SSP585
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';

load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

%%
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';

load([fpath 'CanESM_historical_biomes.mat'],'hcbiomes');
load([fpath 'CNRM_historical_biomes.mat'],'hnbiomes');
load([fpath 'GFDL_historical_biomes.mat'],'hgbiomes');
load([fpath 'IPSL_historical_biomes.mat'],'hibiomes');
load([fpath 'UKESM_historical_biomes.mat'],'hubiomes');

load([fpath 'CanESM_ssp585_biomes.mat'],'scbiomes');
load([fpath 'CNRM_ssp585_biomes.mat'],'snbiomes');
load([fpath 'GFDL_ssp585_biomes.mat'],'sgbiomes');
load([fpath 'IPSL_ssp585_biomes.mat'],'sibiomes');
load([fpath 'UKESM_ssp585_biomes.mat'],'lat','lon','subiomes');

[lat_g,lon_g] = meshgrid(lat,lon);
clear lat lon

%%
load([bpath 'data_biomes_MODISAqua_x1.mat']);
[lat_c,lon_c] = meshgrid(lat,lon);

%% Maps

clatlim=[-90 90];
clonlim=[-280 80];
load coastlines;

cmv = colormap(viridis(3));

%% Get rid of seam
lon_b = lon_g;
lat_b = lat_g;
hcbiomes_s = hcbiomes;
hnbiomes_s = hnbiomes;
hgbiomes_s = hgbiomes;
hibiomes_s = hibiomes;
hubiomes_s = hubiomes;
scbiomes_s = scbiomes;
snbiomes_s = snbiomes;
sgbiomes_s = sgbiomes;
sibiomes_s = sibiomes;
subiomes_s = subiomes;

lon_b(361,:) = lon_b(360,:)+1;
lat_b(361,:) = lat_b(360,:);
hcbiomes_s(361,:) = hcbiomes(360,:);
hnbiomes_s(361,:) = hnbiomes(360,:);
hgbiomes_s(361,:) = hgbiomes(360,:);
hibiomes_s(361,:) = hibiomes(360,:);
hubiomes_s(361,:) = hubiomes(360,:);
scbiomes_s(361,:) = scbiomes(360,:);
snbiomes_s(361,:) = snbiomes(360,:);
sgbiomes_s(361,:) = sgbiomes(360,:);
sibiomes_s(361,:) = sibiomes(360,:);
subiomes_s(361,:) = subiomes(360,:);

%% load troph amp data
load('troph_amp_ssp585_chl_mesoz200_10yr.mat')

%% check same orientation
figure
pcolor(lat_b)
shading flat
title('latb')
figure
pcolor(lat_s)
shading flat
title('latb')
figure
pcolor(lon_b)
shading flat
title('lats')
figure
pcolor(lon_s)
shading flat
title('lons')
figure
pcolor(hcbiomes_s)
shading flat
title('hcbiomes')
figure
pcolor(cta1)
shading flat
title('cta1')

%% flip biomes
hcbiomes_s = fliplr(hcbiomes_s);
hnbiomes_s = fliplr(hnbiomes_s);
hgbiomes_s = fliplr(hgbiomes_s);
hibiomes_s = fliplr(hibiomes_s);
hubiomes_s = fliplr(hubiomes_s);

scbiomes_s = fliplr(scbiomes_s);
snbiomes_s = fliplr(snbiomes_s);
sgbiomes_s = fliplr(sgbiomes_s);
sibiomes_s = fliplr(sibiomes_s);
subiomes_s = fliplr(subiomes_s);

%% Take means by biome
% Exclude LC > 45N/S
pol = find(lat_b <= 45 & lat_b >= -45);

clc = find(hcbiomes_s==1);
nlc = find(hnbiomes_s==1);
glc = find(hgbiomes_s==1);
ilc = find(hibiomes_s==1);
ulc = find(hubiomes_s==1);

cid = intersect(clc,pol);
nid = intersect(nlc,pol);
gid = intersect(glc,pol);
iid = intersect(ilc,pol);
uid = intersect(ulc,pol);

css = (hcbiomes_s==2);
nss = (hnbiomes_s==2);
gss = (hgbiomes_s==2);
iss = (hibiomes_s==2);
uss = (hubiomes_s==2);

cps = (hcbiomes_s==3);
nps = (hnbiomes_s==3);
gps = (hgbiomes_s==3);
ips = (hibiomes_s==3);
ups = (hubiomes_s==3);

%% ssp585
clc2 = find(scbiomes_s==1);
nlc2 = find(snbiomes_s==1);
glc2 = find(sgbiomes_s==1);
ilc2 = find(sibiomes_s==1);
ulc2 = find(subiomes_s==1);

cid2 = intersect(clc2,pol);
nid2 = intersect(nlc2,pol);
gid2 = intersect(glc2,pol);
iid2 = intersect(ilc2,pol);
uid2 = intersect(ulc2,pol);

css2 = (scbiomes_s==2);
nss2 = (snbiomes_s==2);
gss2 = (sgbiomes_s==2);
iss2 = (sibiomes_s==2);
uss2 = (subiomes_s==2);

cps2 = (scbiomes_s==3);
nps2 = (snbiomes_s==3);
gps2 = (sgbiomes_s==3);
ips2 = (sibiomes_s==3);
ups2 = (subiomes_s==3);

%% throw out high outliers
gh=quantile(gta1(:),0.99);
ch=quantile(cta1(:),0.99);

gta1(gta1 > gh) = nan;
cta1(cta1 > ch) = nan;

%% Using hist biomes
mta(1,1) = nanmean(cta1(:));
mta(2,1) = nanmean(nta1(:));
mta(3,1) = nanmean(gta1(:));
mta(4,1) = nanmean(ita1(:));
mta(5,1) = nanmean(uta1(:));

mta(1,2) = nanmean(cta1(cid));
mta(2,2) = nanmean(nta1(nid));
mta(3,2) = nanmean(gta1(gid));
mta(4,2) = nanmean(ita1(iid));
mta(5,2) = nanmean(uta1(uid));

mta(1,3) = nanmean(cta1(css));
mta(2,3) = nanmean(nta1(nss));
mta(3,3) = nanmean(gta1(gss));
mta(4,3) = nanmean(ita1(iss));
mta(5,3) = nanmean(uta1(uss));

mta(1,4) = nanmean(cta1(cps));
mta(2,4) = nanmean(nta1(nps));
mta(3,4) = nanmean(gta1(gps));
mta(4,4) = nanmean(ita1(ips));
mta(5,4) = nanmean(uta1(ups));

%% Using ssp biomes
mta(1,5) = nanmean(cta1(cid2));
mta(2,5) = nanmean(nta1(nid2));
mta(3,5) = nanmean(gta1(gid2));
mta(4,5) = nanmean(ita1(iid2));
mta(5,5) = nanmean(uta1(uid2));

mta(1,6) = nanmean(cta1(css2));
mta(2,6) = nanmean(nta1(nss2));
mta(3,6) = nanmean(gta1(gss2));
mta(4,6) = nanmean(ita1(iss2));
mta(5,6) = nanmean(uta1(uss2));

mta(1,7) = nanmean(cta1(cps2));
mta(2,7) = nanmean(nta1(nps2));
mta(3,7) = nanmean(gta1(gps2));
mta(4,7) = nanmean(ita1(ips2));
mta(5,7) = nanmean(uta1(ups2));

%% table 

Tta = array2table(mta,'RowNames',{'CAN','CNRM','GFDL','IPSL','UK'},...
    'VariableNames',{'Global','hLC','hHCSS','hHCPS','sLC','sHCSS','sHCPS'});

save('troph_amp_ssp585_chl_mesoz200_10yr.mat','mta','Tta','-append');
writetable(Tta,'mean_troph_amp_ssp585_chl_mesoz200_10yr_biomes.csv',...
    'Delimiter',',','WriteRowNames',true)


