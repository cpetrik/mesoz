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
load([fpath 'CMCC_historical_biomes.mat'],'hmbiomes');
load([fpath 'CNRM_historical_biomes.mat'],'hnbiomes');
load([fpath 'GFDL_historical_biomes.mat'],'hgbiomes');
load([fpath 'IPSL_historical_biomes.mat'],'hibiomes');
load([fpath 'UKESM_historical_biomes.mat'],'hubiomes');

load([fpath 'CanESM_ssp585_biomes.mat'],'scbiomes');
load([fpath 'CMCC_ssp585_biomes.mat'],'smbiomes');
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
close all

lon_b = lon_g;
lat_b = lat_g;

%% load troph amp data
load('troph_amp_ssp585_chl_mesoz200_10yr.mat')

% use regular grid
cta1 = cta1(1:360,:);
mta1 = mta1(1:360,:);
nta1 = nta1(1:360,:);
gta1 = gta1(1:360,:);
ita1 = ita1(1:360,:);
uta1 = uta1(1:360,:);

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
title('lonb')
figure
pcolor(lon_s)
shading flat
title('lons')
figure
pcolor(hcbiomes)
shading flat
title('hcbiomes')
figure
pcolor(cta1)
shading flat
title('cta1')

%% flip biomes
close all
hcbiomes = fliplr(hcbiomes);
hmbiomes = fliplr(hmbiomes);
hnbiomes = fliplr(hnbiomes);
hgbiomes = fliplr(hgbiomes);
hibiomes = fliplr(hibiomes);
hubiomes = fliplr(hubiomes);

scbiomes = fliplr(scbiomes);
smbiomes = fliplr(smbiomes);
snbiomes = fliplr(snbiomes);
sgbiomes = fliplr(sgbiomes);
sibiomes = fliplr(sibiomes);
subiomes = fliplr(subiomes);

%% Take means by biome
% Exclude LC > 45N/S
pol = find(lat_b <= 45 & lat_b >= -45);

clc = find(hcbiomes==1);
mlc = find(hmbiomes==1);
nlc = find(hnbiomes==1);
glc = find(hgbiomes==1);
ilc = find(hibiomes==1);
ulc = find(hubiomes==1);

cid = intersect(clc,pol);
mid = intersect(mlc,pol);
nid = intersect(nlc,pol);
gid = intersect(glc,pol);
iid = intersect(ilc,pol);
uid = intersect(ulc,pol);

css = (hcbiomes==2);
mss = (hmbiomes==2);
nss = (hnbiomes==2);
gss = (hgbiomes==2);
iss = (hibiomes==2);
uss = (hubiomes==2);

cps = (hcbiomes==3);
mps = (hmbiomes==3);
nps = (hnbiomes==3);
gps = (hgbiomes==3);
ips = (hibiomes==3);
ups = (hubiomes==3);

%% ssp585
clc2 = find(scbiomes==1);
mlc2 = find(smbiomes==1);
nlc2 = find(snbiomes==1);
glc2 = find(sgbiomes==1);
ilc2 = find(sibiomes==1);
ulc2 = find(subiomes==1);

cid2 = intersect(clc2,pol);
mid2 = intersect(mlc2,pol);
nid2 = intersect(nlc2,pol);
gid2 = intersect(glc2,pol);
iid2 = intersect(ilc2,pol);
uid2 = intersect(ulc2,pol);

css2 = (scbiomes==2);
mss2 = (smbiomes==2);
nss2 = (snbiomes==2);
gss2 = (sgbiomes==2);
iss2 = (sibiomes==2);
uss2 = (subiomes==2);

cps2 = (scbiomes==3);
mps2 = (smbiomes==3);
nps2 = (snbiomes==3);
gps2 = (sgbiomes==3);
ips2 = (sibiomes==3);
ups2 = (subiomes==3);

%% throw out high outliers
gh=quantile(gta1(:),0.99);
ch=quantile(cta1(:),0.99);

gta1(gta1 > gh) = nan;
cta1(cta1 > ch) = nan;

%% Using hist biomes
mta(1,1) = nansum(cta1(:) .* cell_area(:)) ./ nansum(cell_area(:));
mta(2,1) = nansum(mta1(:) .* cell_area(:)) ./ nansum(cell_area(:));
mta(3,1) = nansum(nta1(:) .* cell_area(:)) ./ nansum(cell_area(:));
mta(4,1) = nansum(gta1(:) .* cell_area(:)) ./ nansum(cell_area(:));
mta(5,1) = nansum(ita1(:) .* cell_area(:)) ./ nansum(cell_area(:));
mta(6,1) = nansum(uta1(:) .* cell_area(:)) ./ nansum(cell_area(:));

mta(1,2) = nansum(cta1(cid) .* cell_area(cid)) ./ nansum(cell_area(cid));
mta(2,2) = nansum(mta1(mid) .* cell_area(mid)) ./ nansum(cell_area(mid));
mta(3,2) = nansum(nta1(nid) .* cell_area(nid)) ./ nansum(cell_area(nid));
mta(4,2) = nansum(gta1(gid) .* cell_area(gid)) ./ nansum(cell_area(gid));
mta(5,2) = nansum(ita1(iid) .* cell_area(iid)) ./ nansum(cell_area(iid));
mta(6,2) = nansum(uta1(uid) .* cell_area(uid)) ./ nansum(cell_area(uid));

mta(1,3) = nansum(cta1(css) .* cell_area(css)) ./ nansum(cell_area(css));
mta(2,3) = nansum(mta1(mss) .* cell_area(mss)) ./ nansum(cell_area(mss));
mta(3,3) = nansum(nta1(nss) .* cell_area(nss)) ./ nansum(cell_area(nss));
mta(4,3) = nansum(gta1(gss) .* cell_area(gss)) ./ nansum(cell_area(gss));
mta(5,3) = nansum(ita1(iss) .* cell_area(iss)) ./ nansum(cell_area(iss));
mta(6,3) = nansum(uta1(uss) .* cell_area(uss)) ./ nansum(cell_area(uss));

mta(1,4) = nansum(cta1(cps) .* cell_area(cps)) ./ nansum(cell_area(cps));
mta(2,4) = nansum(mta1(mps) .* cell_area(mps)) ./ nansum(cell_area(mps));
mta(3,4) = nansum(nta1(nps) .* cell_area(nps)) ./ nansum(cell_area(nps));
mta(4,4) = nansum(gta1(gps) .* cell_area(gps)) ./ nansum(cell_area(gps));
mta(5,4) = nansum(ita1(ips) .* cell_area(ips)) ./ nansum(cell_area(ips));
mta(6,4) = nansum(uta1(ups) .* cell_area(ups)) ./ nansum(cell_area(ups));

%% Using ssp biomes
mta(1,5) = nansum(cta1(cid2) .* cell_area(cid2)) ./ nansum(cell_area(cid2));
mta(2,5) = nansum(mta1(mid2) .* cell_area(mid2)) ./ nansum(cell_area(mid2));
mta(3,5) = nansum(nta1(nid2) .* cell_area(nid2)) ./ nansum(cell_area(nid2));
mta(4,5) = nansum(gta1(gid2) .* cell_area(gid2)) ./ nansum(cell_area(gid2));
mta(5,5) = nansum(ita1(iid2) .* cell_area(iid2)) ./ nansum(cell_area(iid2));
mta(6,5) = nansum(uta1(uid2) .* cell_area(uid2)) ./ nansum(cell_area(uid2));

mta(1,6) = nansum(cta1(css2) .* cell_area(css2)) ./ nansum(cell_area(css2));
mta(2,6) = nansum(mta1(mss2) .* cell_area(mss2)) ./ nansum(cell_area(mss2));
mta(3,6) = nansum(nta1(nss2) .* cell_area(nss2)) ./ nansum(cell_area(nss2));
mta(4,6) = nansum(gta1(gss2) .* cell_area(gss2)) ./ nansum(cell_area(gss2));
mta(5,6) = nansum(ita1(iss2) .* cell_area(iss2)) ./ nansum(cell_area(iss2));
mta(6,6) = nansum(uta1(uss2) .* cell_area(uss2)) ./ nansum(cell_area(uss2));

mta(1,7) = nansum(cta1(cps2) .* cell_area(cps2)) ./ nansum(cell_area(cps2));
mta(2,7) = nansum(mta1(mps2) .* cell_area(mps2)) ./ nansum(cell_area(mps2));
mta(3,7) = nansum(nta1(nps2) .* cell_area(nps2)) ./ nansum(cell_area(nps2));
mta(4,7) = nansum(gta1(gps2) .* cell_area(gps2)) ./ nansum(cell_area(gps2));
mta(5,7) = nansum(ita1(ips2) .* cell_area(ips2)) ./ nansum(cell_area(ips2));
mta(6,7) = nansum(uta1(ups2) .* cell_area(ups2)) ./ nansum(cell_area(ups2));

%% table 

Tta = array2table(mta,'RowNames',{'CAN','CMCC','CNRM','GFDL','IPSL','UK'},...
    'VariableNames',{'Global','hLC','hHCSS','hHCPS','sLC','sHCSS','sHCPS'});

save('troph_amp_ssp585_chl_mesoz200_10yr_areaw.mat','mta','Tta');
writetable(Tta,'mean_troph_amp_ssp585_chl_mesoz200_10yr_biomes_areaw.csv',...
    'Delimiter',',','WriteRowNames',true)


