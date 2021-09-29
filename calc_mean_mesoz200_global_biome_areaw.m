% CMIP6 output 
% 200m integrations
% Area-weighted means globally and by biome

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';

load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

%% load fixed
load('Hist_mesoz_clim_biomes_samegrid.mat')

clear glmz_mo

%% GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load('glm100_obs_mesoz.mat');
load('glm100_obs_grid.mat');

% Reshape
[ni,nj,nt] = size(gz);
oz_mo = reshape(glmobsmesoz,ni,nj,12);
lat_o = reshape(Lat,ni,nj);
lon_o = reshape(Lon,ni,nj);

%% Check new glmm orientation
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

%% Reshape
close all

cz = reshape(cz,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
nz = reshape(nz,ni*nj,nt);
gz = reshape(gz,ni*nj,nt);
iz = reshape(iz,ni*nj,nt);
uz = reshape(uz,ni*nj,nt);
oz = reshape(oz_mo,ni*nj,nt);

area = repmat(cell_area,1,1,nt);
areav = reshape(area,ni*nj,nt);

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

%% Isolate
Czlc = (cz(cid,:) .* areav(cid,:));
Mzlc = (mz(mid,:) .* areav(mid,:));
Nzlc = (nz(nid,:) .* areav(nid,:));
Gzlc = (gz(gid,:) .* areav(gid,:));
Izlc = (iz(iid,:) .* areav(iid,:));
Uzlc = (uz(uid,:) .* areav(uid,:));
Ozlc = (oz(oid,:) .* areav(oid,:));

Calc = (areav(cid,:));
Malc = (areav(mid,:));
Nalc = (areav(nid,:));
Galc = (areav(gid,:));
Ialc = (areav(iid,:));
Ualc = (areav(uid,:));
Oalc = (areav(oid,:));

Czss = (cz(css,:) .* areav(css,:));
Mzss = (mz(mss,:) .* areav(mss,:));
Nzss = (nz(nss,:) .* areav(nss,:));
Gzss = (gz(gss,:) .* areav(gss,:));
Izss = (iz(iss,:) .* areav(iss,:));
Uzss = (uz(uss,:) .* areav(uss,:));
Ozss = (oz(oss,:) .* areav(oss,:));

Cass = (areav(css,:));
Mass = (areav(mss,:));
Nass = (areav(nss,:));
Gass = (areav(gss,:));
Iass = (areav(iss,:));
Uass = (areav(uss,:));
Oass = (areav(oss,:));

Czps = (cz(cps,:) .* areav(cps,:));
Mzps = (mz(mps,:) .* areav(mps,:));
Nzps = (nz(nps,:) .* areav(nps,:));
Gzps = (gz(gps,:) .* areav(gps,:));
Izps = (iz(ips,:) .* areav(ips,:));
Uzps = (uz(ups,:) .* areav(ups,:));
Ozps = (oz(ops,:) .* areav(ops,:));

Caps = (areav(cps,:));
Maps = (areav(mps,:));
Naps = (areav(nps,:));
Gaps = (areav(gps,:));
Iaps = (areav(ips,:));
Uaps = (areav(ups,:));
Oaps = (areav(ops,:));

%% Take means by biome
zmeans(1,1) = nansum(cz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(2,1) = nansum(mz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(3,1) = nansum(nz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(4,1) = nansum(gz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(5,1) = nansum(iz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(6,1) = nansum(uz(:) .* areav(:)) ./ nansum(areav(:));
zmeans(7,1) = nansum(oz(:) .* areav(:)) ./ nansum(areav(:));

% LC
zmeans(1,2) = nansum(Czlc(:)) ./ nansum(Calc(:));
zmeans(2,2) = nansum(Mzlc(:)) ./ nansum(Malc(:));
zmeans(3,2) = nansum(Nzlc(:)) ./ nansum(Nalc(:));
zmeans(4,2) = nansum(Gzlc(:)) ./ nansum(Galc(:));
zmeans(5,2) = nansum(Izlc(:)) ./ nansum(Ialc(:));
zmeans(6,2) = nansum(Uzlc(:)) ./ nansum(Ualc(:));
zmeans(7,2) = nansum(Ozlc(:)) ./ nansum(Oalc(:));
% SS
zmeans(1,3) = nansum(Czss(:)) ./ nansum(Cass(:));
zmeans(2,3) = nansum(Mzss(:)) ./ nansum(Mass(:));
zmeans(3,3) = nansum(Nzss(:)) ./ nansum(Nass(:));
zmeans(4,3) = nansum(Gzss(:)) ./ nansum(Gass(:));
zmeans(5,3) = nansum(Izss(:)) ./ nansum(Iass(:));
zmeans(6,3) = nansum(Uzss(:)) ./ nansum(Uass(:));
zmeans(7,3) = nansum(Ozss(:)) ./ nansum(Oass(:));
% PS
zmeans(1,4) = nansum(Czps(:)) ./ nansum(Caps(:));
zmeans(2,4) = nansum(Mzps(:)) ./ nansum(Maps(:));
zmeans(3,4) = nansum(Nzps(:)) ./ nansum(Naps(:));
zmeans(4,4) = nansum(Gzps(:)) ./ nansum(Gaps(:));
zmeans(5,4) = nansum(Izps(:)) ./ nansum(Iaps(:));
zmeans(6,4) = nansum(Uzps(:)) ./ nansum(Uaps(:));
zmeans(7,4) = nansum(Ozps(:)) ./ nansum(Oaps(:));

%%
simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obs'};
sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

Tmeanz = array2table(zmeans,'RowNames',simtext,'VariableNames',...
    {'Global','LC','HCSS','HCPS'});
writetable(Tmeanz,[sfile 'means_areaw_hist_aclim_zmeso200_obsglm100_global_biomes.csv'],'WriteRowNames',true);
save('means_areaw_hist_aclim_zmeso200_obsglm100_global_biomes.mat',...
    'Tmeanz','zmeans');












