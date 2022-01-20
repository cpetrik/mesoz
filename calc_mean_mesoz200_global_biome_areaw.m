% CMIP6 output 
% 200m integrations
% Area-weighted means globally and by biome

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';

load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

%% load fixed ESMs
load('Hist_mesoz_clim_biomes_samegrid.mat')

clear glmz_mo

%% GLM zmeso, grid vars
% opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
% load('glm100_obs_mesoz.mat');
% load('glm100_obs_grid.mat');
% 
% % Reshape
[ni,nj,nt] = size(gz);
% oz_mo = reshape(glmobsmesoz,ni,nj,12);
% lat_o = reshape(Lat,ni,nj);
% lon_o = reshape(Lon,ni,nj);
% 
% %% Stromberg obs model
% fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/Stromberg_x1_all/';
% load([fpath 'StrombergQTR_clim_int200_mgCm2.mat']);
% [lat_m,lon_m] = meshgrid(lat,lon);
% 
% sz_mo = sz;
% 
% clear lat lon sz
% 
% Biomes
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';

load([bpath 'SeaWiFS_based_biomes_x1.mat']);
[lat_b,lon_b] = meshgrid(lat,lon);

clear lat lon
% 
% %% COPEPOD
% load('copepod-2012_cmass_all_gridded.mat','lat_g','lon_g',...
%     'zoo_g','fileid','units')
% 
% % From m-3 to m-2 (integrate top 200m)
% zoo_200 = zoo_g*200;

%% check orientations
figure(3)
pcolor(lat_m); shading flat;
title('Zlat')

figure(4)
pcolor(lat_b); shading flat;
title('Blat')

%% Check new glmm orientation
figure(1)
pcolor(squeeze(sz_mo(:,:,6))); shading flat;
title('Obs SM')

figure(5)
pcolor(zoo_200); shading flat;
title('Obs C')

figure(2)
pcolor(sbiomes); shading flat;
title('obsSM')

figure
pcolor(gz(:,:,6)); shading flat;
title('GFDL')

figure
pcolor(oz_mo(:,:,6)); shading flat;
title('obsGLMM')

figure
pcolor(hgbiomes); shading flat;
title('G biome')

figure
pcolor(biomes); shading flat;
title('obs biome')

%% Flip
% close all
% 
% sz_mo = fliplr(sz_mo);
% lat_m = fliplr(lat_m);
% 
% save('Hist_mesoz_clim_biomes_samegrid.mat','oz_mo','sz_mo','zoo_200','sbiomes',...
%     '-append')

%% Reshape
close all

cz = reshape(cz,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
nz = reshape(nz,ni*nj,nt);
gz = reshape(gz,ni*nj,nt);
iz = reshape(iz,ni*nj,nt);
uz = reshape(uz,ni*nj,nt);
oz = reshape(oz_mo,ni*nj,nt);
sz = reshape(sz_mo,ni*nj,nt);
z2 = reshape(zoo_200,ni*nj,1);

area = repmat(cell_area,1,1,nt);
areav = reshape(area,ni*nj,nt);

znan = isnan(zoo_200);
area2 = cell_area;
area2(znan) = nan;
areav2 = reshape(area2,ni*nj,1);

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
slc = find(sbiomes==1);

cid = intersect(clc,pol);
mid = intersect(mlc,pol);
nid = intersect(nlc,pol);
gid = intersect(glc,pol);
iid = intersect(ilc,pol);
uid = intersect(ulc,pol);
oid = intersect(olc,pol);
sid = intersect(slc,pol);

css = (hcbiomes==2);
mss = (hmbiomes==2);
nss = (hnbiomes==2);
gss = (hgbiomes==2);
iss = (hibiomes==2);
uss = (hubiomes==2);
oss = (biomes==2);
sss = (sbiomes==2);

cps = (hcbiomes==3);
mps = (hmbiomes==3);
nps = (hnbiomes==3);
gps = (hgbiomes==3);
ips = (hibiomes==3);
ups = (hubiomes==3);
ops = (biomes==3);
sps = (sbiomes==3);

%% Isolate
Czlc = (cz(cid,:) .* areav(cid,:));
Mzlc = (mz(mid,:) .* areav(mid,:));
Nzlc = (nz(nid,:) .* areav(nid,:));
Gzlc = (gz(gid,:) .* areav(gid,:));
Izlc = (iz(iid,:) .* areav(iid,:));
Uzlc = (uz(uid,:) .* areav(uid,:));
Ozlc = (oz(oid,:) .* areav(oid,:));
Szlc = (sz(sid,:) .* areav(sid,:));
CMzlc = (z2(oid,:) .* areav2(oid,:));
CSzlc = (z2(sid,:) .* areav2(sid,:));

Calc = (areav(cid,:));
Malc = (areav(mid,:));
Nalc = (areav(nid,:));
Galc = (areav(gid,:));
Ialc = (areav(iid,:));
Ualc = (areav(uid,:));
Oalc = (areav(oid,:));
Salc = (areav(sid,:));
CMalc = (areav2(oid,:));
CSalc = (areav2(sid,:));

Czss = (cz(css,:) .* areav(css,:));
Mzss = (mz(mss,:) .* areav(mss,:));
Nzss = (nz(nss,:) .* areav(nss,:));
Gzss = (gz(gss,:) .* areav(gss,:));
Izss = (iz(iss,:) .* areav(iss,:));
Uzss = (uz(uss,:) .* areav(uss,:));
Ozss = (oz(oss,:) .* areav(oss,:));
Szss = (sz(sss,:) .* areav(sss,:));
CMzss = (z2(oss,:) .* areav2(oss,:));
CSzss = (z2(sss,:) .* areav2(sss,:));

Cass = (areav(css,:));
Mass = (areav(mss,:));
Nass = (areav(nss,:));
Gass = (areav(gss,:));
Iass = (areav(iss,:));
Uass = (areav(uss,:));
Oass = (areav(oss,:));
Sass = (areav(sss,:));
CMass = (areav2(oss,:));
CSass = (areav2(sss,:));

Czps = (cz(cps,:) .* areav(cps,:));
Mzps = (mz(mps,:) .* areav(mps,:));
Nzps = (nz(nps,:) .* areav(nps,:));
Gzps = (gz(gps,:) .* areav(gps,:));
Izps = (iz(ips,:) .* areav(ips,:));
Uzps = (uz(ups,:) .* areav(ups,:));
Ozps = (oz(ops,:) .* areav(ops,:));
Szps = (sz(sps,:) .* areav(sps,:));
CMzps = (z2(ops,:) .* areav2(ops,:));
CSzps = (z2(sps,:) .* areav2(sps,:));

Caps = (areav(cps,:));
Maps = (areav(mps,:));
Naps = (areav(nps,:));
Gaps = (areav(gps,:));
Iaps = (areav(ips,:));
Uaps = (areav(ups,:));
Oaps = (areav(ops,:));
Saps = (areav(sps,:));
CMaps = (areav2(ops,:));
CSaps = (areav2(sps,:));

%% Take means by biome
zmeans(1,1) = nansum(cz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(2,1) = nansum(mz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(3,1) = nansum(nz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(4,1) = nansum(gz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(5,1) = nansum(iz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(6,1) = nansum(uz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(7,1) = nansum(oz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(8,1) = nansum(sz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(9,1) = nansum(z2(:) .* areav2(:)) ./ nansum(areav2(:));
zmeans(10,1) = nansum(z2(:) .* areav2(:)) ./ nansum(areav2(:));

% LC
zmeans(1,2) = nansum(Czlc(:)) ./ nansum(Calc(:));
zmeans(2,2) = nansum(Mzlc(:)) ./ nansum(Malc(:));
zmeans(3,2) = nansum(Nzlc(:)) ./ nansum(Nalc(:));
zmeans(4,2) = nansum(Gzlc(:)) ./ nansum(Galc(:));
zmeans(5,2) = nansum(Izlc(:)) ./ nansum(Ialc(:));
zmeans(6,2) = nansum(Uzlc(:)) ./ nansum(Ualc(:));
zmeans(7,2) = nansum(Ozlc(:)) ./ nansum(Oalc(:));
zmeans(8,2) = nansum(Szlc(:)) ./ nansum(Salc(:));
zmeans(9,2) = nansum(CMzlc(:)) ./ nansum(CMalc(:));
zmeans(10,2) = nansum(CSzlc(:)) ./ nansum(CSalc(:));
% SS
zmeans(1,3) = nansum(Czss(:)) ./ nansum(Cass(:));
zmeans(2,3) = nansum(Mzss(:)) ./ nansum(Mass(:));
zmeans(3,3) = nansum(Nzss(:)) ./ nansum(Nass(:));
zmeans(4,3) = nansum(Gzss(:)) ./ nansum(Gass(:));
zmeans(5,3) = nansum(Izss(:)) ./ nansum(Iass(:));
zmeans(6,3) = nansum(Uzss(:)) ./ nansum(Uass(:));
zmeans(7,3) = nansum(Ozss(:)) ./ nansum(Oass(:));
zmeans(8,3) = nansum(Szss(:)) ./ nansum(Sass(:));
zmeans(9,3) = nansum(CMzss(:)) ./ nansum(CMass(:));
zmeans(10,3) = nansum(CSzss(:)) ./ nansum(CSass(:));
% PS
zmeans(1,4) = nansum(Czps(:)) ./ nansum(Caps(:));
zmeans(2,4) = nansum(Mzps(:)) ./ nansum(Maps(:));
zmeans(3,4) = nansum(Nzps(:)) ./ nansum(Naps(:));
zmeans(4,4) = nansum(Gzps(:)) ./ nansum(Gaps(:));
zmeans(5,4) = nansum(Izps(:)) ./ nansum(Iaps(:));
zmeans(6,4) = nansum(Uzps(:)) ./ nansum(Uaps(:));
zmeans(7,4) = nansum(Ozps(:)) ./ nansum(Oaps(:));
zmeans(8,4) = nansum(Szps(:)) ./ nansum(Saps(:));
zmeans(9,4) = nansum(CMzps(:)) ./ nansum(CMaps(:));
zmeans(10,4) = nansum(CSzps(:)) ./ nansum(CSaps(:));

%%
simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLMM','obsSM','obsCM','obsCS'};
sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

Tmeanz = array2table(zmeans,'RowNames',simtext,'VariableNames',...
    {'Global','LC','HCSS','HCPS'});

%%
writetable(Tmeanz,[sfile 'means_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.csv'],'WriteRowNames',true);
save('means_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.mat',...
    'Tmeanz','zmeans');

%% test COPE not area-weighted
%area-weighting lowering the means b/c not all grid cells with area have
%data
cmeans(1,1) = nanmean(z2(:));
cmeans(2,1) = nanmean(z2(:));
% LC
cmeans(1,2) = nanmean(z2(oid,:));
cmeans(2,2) = nanmean(z2(sid,:));
% SS
cmeans(1,3) = nanmean(z2(oss,:));
cmeans(2,3) = nanmean(z2(sss,:));
% PS
cmeans(1,4) = nanmean(z2(ops,:));
cmeans(2,4) = nanmean(z2(sps,:));








