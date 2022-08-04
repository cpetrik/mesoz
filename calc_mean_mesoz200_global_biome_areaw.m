% CMIP6 output 
% 200m integrations
% Area-weighted means globally and by biome

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
apath='/Volumes/MIP/Fish-MIP/CMIP6/';
%apath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/';

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
%bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/biome_masks/';

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

Tstdz = array2table(zmeans,'RowNames',simtext,'VariableNames',...
    {'Global','LC','HCSS','HCPS'});

%%
writetable(Tstdz,[sfile 'means_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.csv'],'WriteRowNames',true);
save('means_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.mat',...
    'Tmeanz','zmeans');

%% area-weighted std dev
% num = sum (area .* (obs - weightedmean).^2)
% den = ((n-1)*sum(area))/n
% wsd = sqrt( num / den )

% area-weighted means
czwm = zmeans(1,1) ;
mzwm = zmeans(2,1) ;
nzwm = zmeans(3,1) ;
gzwm = zmeans(4,1) ;
izwm = zmeans(5,1) ;
uzwm = zmeans(6,1) ;
ozwm = zmeans(7,1) ;
szwm = zmeans(8,1) ;
cmzwm = zmeans(9,1) ;
cszwm = zmeans(10,1);

% Numerator
czg = nansum(areav(:) .* ((cz(:) - czwm).^2));
mzg = nansum(areav(:) .* ((mz(:) - mzwm).^2));
nzg = nansum(areav(:) .* ((nz(:) - nzwm).^2));
gzg = nansum(areav(:) .* ((gz(:) - gzwm).^2));
izg = nansum(areav(:) .* ((iz(:) - izwm).^2));
uzg = nansum(areav(:) .* ((uz(:) - uzwm).^2));
ozg = nansum(areav(:) .* ((oz(:) - ozwm).^2));
szg = nansum(areav(:) .* ((sz(:) - szwm).^2));
cmzg = nansum(areav2(:) .* ((z2(:) - cmzwm).^2));
cszg = nansum(areav2(:) .* ((z2(:) - cszwm).^2));

%Szlc = (sz(sid,:) .* areav(sid,:));
%CMzlc = (z2(oid,:) .* areav2(oid,:));
Czlc = nansum(areav(cid,:) .* ((cz(cid,:) - czwm).^2));
Mzlc = nansum(areav(mid,:) .* ((mz(mid,:) - mzwm).^2));
Nzlc = nansum(areav(nid,:) .* ((nz(nid,:) - nzwm).^2));
Gzlc = nansum(areav(gid,:) .* ((gz(gid,:) - gzwm).^2));
Izlc = nansum(areav(iid,:) .* ((iz(iid,:) - izwm).^2));
Uzlc = nansum(areav(uid,:) .* ((uz(uid,:) - uzwm).^2));
Ozlc = nansum(areav(oid,:) .* ((oz(oid,:) - ozwm).^2));
Szlc = nansum(areav(sid,:) .* ((sz(sid,:) - szwm).^2));
CMzlc = nansum(areav2(oid,:) .* ((z2(oid,:) - cmzwm).^2));
CSzlc = nansum(areav2(sid,:) .* ((z2(sid,:) - cszwm).^2));

Czss = nansum(areav(css,:) .* ((cz(css,:) - czwm).^2));
Mzss = nansum(areav(mss,:) .* ((mz(mss,:) - mzwm).^2));
Nzss = nansum(areav(nss,:) .* ((nz(nss,:) - nzwm).^2));
Gzss = nansum(areav(gss,:) .* ((gz(gss,:) - gzwm).^2));
Izss = nansum(areav(iss,:) .* ((iz(iss,:) - izwm).^2));
Uzss = nansum(areav(uss,:) .* ((uz(uss,:) - uzwm).^2));
Ozss = nansum(areav(oss,:) .* ((oz(oss,:) - ozwm).^2));
Szss = nansum(areav(sss,:) .* ((sz(sss,:) - szwm).^2));
CMzss = nansum(areav2(oss,:) .* ((z2(oss,:) - cmzwm).^2));
CSzss = nansum(areav2(sss,:) .* ((z2(sss,:) - cszwm).^2));

Czps = nansum(areav(cps,:) .* ((cz(cps,:) - czwm).^2));
Mzps = nansum(areav(mps,:) .* ((mz(mps,:) - mzwm).^2));
Nzps = nansum(areav(nps,:) .* ((nz(nps,:) - nzwm).^2));
Gzps = nansum(areav(gps,:) .* ((gz(gps,:) - gzwm).^2));
Izps = nansum(areav(ips,:) .* ((iz(ips,:) - izwm).^2));
Uzps = nansum(areav(ups,:) .* ((uz(ups,:) - uzwm).^2));
Ozps = nansum(areav(ops,:) .* ((oz(ops,:) - ozwm).^2));
Szps = nansum(areav(sps,:) .* ((sz(sps,:) - szwm).^2));
CMzps = nansum(areav2(ops,:) .* ((z2(ops,:) - cmzwm).^2));
CSzps = nansum(areav2(sps,:) .* ((z2(sps,:) - cszwm).^2));

%% num obs (n)
nclc = 12 * length(cid);
ncss = 12 * sum((css(:)));
ncps = 12 * sum((cps(:)));
nnlc = 12 * length(nid);
nnss = 12 * sum((nss(:)));
nnps = 12 * sum((nps(:)));
nmlc = 12 * length(mid);
nmss = 12 * sum((mss(:)));
nmps = 12 * sum((mps(:)));
nglc = 12 * length(gid);
ngss = 12 * sum((gss(:)));
ngps = 12 * sum((gps(:)));
nilc = 12 * length(iid);
niss = 12 * sum((iss(:)));
nips = 12 * sum((ips(:)));
nulc = 12 * length(uid);
nuss = 12 * sum((uss(:)));
nups = 12 * sum((ups(:)));
nolc = 12 * length(oid);
noss = 12 * sum((oss(:)));
nops = 12 * sum((ops(:)));
nslc = 12 * length(sid);
nsss = 12 * sum((sss(:)));
nsps = 12 * sum((sps(:)));
ncmlc = 12 * length(oid);
ncmss = 12 * sum((oss(:)));
ncmps = 12 * sum((ops(:)));
ncslc = 12 * length(sid);
ncsss = 12 * sum((sss(:)));
ncsps = 12 * sum((sps(:)));

ncz = sum(~isnan(cz(:)));
nnz = sum(~isnan(nz(:)));
nmz = sum(~isnan(mz(:)));
ngz = sum(~isnan(gz(:)));
niz = sum(~isnan(iz(:)));
nuz = sum(~isnan(uz(:)));
noz = sum(~isnan(oz(:)));
nsz = sum(~isnan(sz(:)));
ncmz = sum(~isnan(z2(:)));
ncsz = sum(~isnan(z2(:)));

%% denom = (n-1)*sum(area)/n
cgarea   = sum(areav(:),'omitnan') * (ncz-1) / ncz;
ngarea   = sum(areav(:),'omitnan') * (nnz-1) / nnz;
mgarea   = sum(areav(:),'omitnan') * (nmz-1) / nmz;
ggarea   = sum(areav(:),'omitnan') * (ngz-1) / ngz;
igarea   = sum(areav(:),'omitnan') * (niz-1) / niz;
ugarea   = sum(areav(:),'omitnan') * (nuz-1) / nuz;
ogarea   = sum(areav(:),'omitnan') * (noz-1) / noz;
sgarea   = sum(areav(:),'omitnan') * (nsz-1) / nsz;
cmgarea = sum(areav2(:),'omitnan') * (ncmz-1) / ncmz;
csgarea = sum(areav2(:),'omitnan') * (ncsz-1) / ncsz;

Calc  = sum(areav(cid,:),'omitnan') * (nclc-1) / nclc;
Malc  = sum(areav(mid,:),'omitnan') * (nmlc-1) / nmlc;
Nalc  = sum(areav(nid,:),'omitnan') * (nnlc-1) / nnlc;
Galc  = sum(areav(gid,:),'omitnan') * (nglc-1) / nglc;
Ialc  = sum(areav(iid,:),'omitnan') * (nilc-1) / nilc;
Ualc  = sum(areav(uid,:),'omitnan') * (nulc-1) / nulc;
Oalc  = sum(areav(oid,:),'omitnan') * (nolc-1) / nolc;
Salc  = sum(areav(sid,:),'omitnan') * (nslc-1) / nslc;
CMalc = sum(areav2(oid,:),'omitnan') * (ncmlc-1) / ncmlc;
CSalc = sum(areav2(sid,:),'omitnan') * (ncslc-1) / ncslc;

Cass  = sum(areav(css,:),'omitnan') * (ncss-1) / ncss;
Mass  = sum(areav(mss,:),'omitnan') * (nmss-1) / nmss;
Nass  = sum(areav(nss,:),'omitnan') * (nnss-1) / nnss;
Gass  = sum(areav(gss,:),'omitnan') * (ngss-1) / ngss;
Iass  = sum(areav(iss,:),'omitnan') * (niss-1) / niss;
Uass  = sum(areav(uss,:),'omitnan') * (nuss-1) / nuss;
Oass  = sum(areav(oss,:),'omitnan') * (noss-1) / noss;
Sass  = sum(areav(sss,:),'omitnan') * (nsss-1) / nsss;
CMass = sum(areav2(oss,:),'omitnan') * (ncmss-1) / ncmss;
CSass = sum(areav2(sss,:),'omitnan') * (ncsss-1) / ncsss;

Caps  = sum(areav(cps,:),'omitnan') * (ncps-1) / ncps;
Maps  = sum(areav(mps,:),'omitnan') * (nmps-1) / nmps;
Naps  = sum(areav(nps,:),'omitnan') * (nnps-1) / nnps;
Gaps  = sum(areav(gps,:),'omitnan') * (ngps-1) / ngps;
Iaps  = sum(areav(ips,:),'omitnan') * (nips-1) / nips;
Uaps  = sum(areav(ups,:),'omitnan') * (nups-1) / nups;
Oaps  = sum(areav(ops,:),'omitnan') * (nops-1) / nops;
Saps  = sum(areav(sps,:),'omitnan') * (nsps-1) / nsps;
CMaps = sum(areav2(ops,:),'omitnan') * (ncmps-1) / ncmps;
CSaps = sum(areav2(sps,:),'omitnan') * (ncsps-1) / ncsps;


%% Take s.d. by biome
%global
zstdev(1,1) = sqrt(czg ./ cgarea);
zstdev(2,1) = sqrt(mzg ./ mgarea);
zstdev(3,1) = sqrt(nzg ./ ngarea);
zstdev(4,1) = sqrt(gzg ./ ggarea);
zstdev(5,1) = sqrt(izg ./ igarea);
zstdev(6,1) = sqrt(uzg ./ ugarea);
zstdev(7,1) = sqrt(ozg ./ ogarea);
zstdev(8,1) = sqrt(szg ./ sgarea);
zstdev(9,1) = sqrt(cmzg ./ cmgarea);
zstdev(10,1) = sqrt(cszg ./ csgarea);

% LC
zstdev(1,2) = sqrt(nansum(Czlc(:)) ./ nansum(Calc(:)));
zstdev(2,2) = sqrt(nansum(Mzlc(:)) ./ nansum(Malc(:)));
zstdev(3,2) = sqrt(nansum(Nzlc(:)) ./ nansum(Nalc(:)));
zstdev(4,2) = sqrt(nansum(Gzlc(:)) ./ nansum(Galc(:)));
zstdev(5,2) = sqrt(nansum(Izlc(:)) ./ nansum(Ialc(:)));
zstdev(6,2) = sqrt(nansum(Uzlc(:)) ./ nansum(Ualc(:)));
zstdev(7,2) = sqrt(nansum(Ozlc(:)) ./ nansum(Oalc(:)));
zstdev(8,2) = sqrt(nansum(Szlc(:)) ./ nansum(Salc(:)));
zstdev(9,2) = sqrt(nansum(CMzlc(:)) ./ nansum(CMalc(:)));
zstdev(10,2) = sqrt(nansum(CSzlc(:)) ./ nansum(CSalc(:)));
% SS
zstdev(1,3) = sqrt(nansum(Czss(:)) ./ nansum(Cass(:)));
zstdev(2,3) = sqrt(nansum(Mzss(:)) ./ nansum(Mass(:)));
zstdev(3,3) = sqrt(nansum(Nzss(:)) ./ nansum(Nass(:)));
zstdev(4,3) = sqrt(nansum(Gzss(:)) ./ nansum(Gass(:)));
zstdev(5,3) = sqrt(nansum(Izss(:)) ./ nansum(Iass(:)));
zstdev(6,3) = sqrt(nansum(Uzss(:)) ./ nansum(Uass(:)));
zstdev(7,3) = sqrt(nansum(Ozss(:)) ./ nansum(Oass(:)));
zstdev(8,3) = sqrt(nansum(Szss(:)) ./ nansum(Sass(:)));
zstdev(9,3) = sqrt(nansum(CMzss(:)) ./ nansum(CMass(:)));
zstdev(10,3) = sqrt(nansum(CSzss(:)) ./ nansum(CSass(:)));
% PS
zstdev(1,4) = sqrt(nansum(Czps(:)) ./ nansum(Caps(:)));
zstdev(2,4) = sqrt(nansum(Mzps(:)) ./ nansum(Maps(:)));
zstdev(3,4) = sqrt(nansum(Nzps(:)) ./ nansum(Naps(:)));
zstdev(4,4) = sqrt(nansum(Gzps(:)) ./ nansum(Gaps(:)));
zstdev(5,4) = sqrt(nansum(Izps(:)) ./ nansum(Iaps(:)));
zstdev(6,4) = sqrt(nansum(Uzps(:)) ./ nansum(Uaps(:)));
zstdev(7,4) = sqrt(nansum(Ozps(:)) ./ nansum(Oaps(:)));
zstdev(8,4) = sqrt(nansum(Szps(:)) ./ nansum(Saps(:)));
zstdev(9,4) = sqrt(nansum(CMzps(:)) ./ nansum(CMaps(:)));
zstdev(10,4) = sqrt(nansum(CSzps(:)) ./ nansum(CSaps(:)));

%%
Tstdz = array2table(zstdev,'RowNames',simtext,'VariableNames',...
    {'Global','LC','HCSS','HCPS'});

%%
writetable(Tstdz,[sfile 'stds_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.csv'],'WriteRowNames',true);
save('stds_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.mat',...
    'Tstdz','zstdev');

%% rought test
mean(gz(:),'omitnan')
nanmean(nanmean(gz(gid,:)))
nanmean(nanmean(gz(gss,:)))
nanmean(nanmean(gz(gps,:)))

gtestLC = gz(gid,:);
gtestSS = gz(gss,:);
gtestPS = gz(gps,:);
std(gz(:),'omitnan')
std(gtestLC(:),'omitnan')
std(gtestSS(:),'omitnan')
std(gtestPS(:),'omitnan')

%% Put means and stddevs in one table
tboth(:,1) = zmeans(:,1);
tboth(:,2) = zstdev(:,1);
tboth(:,3) = zmeans(:,2);
tboth(:,4) = zstdev(:,2);
tboth(:,5) = zmeans(:,3);
tboth(:,6) = zstdev(:,3);
tboth(:,7) = zmeans(:,4);
tboth(:,8) = zstdev(:,4);

Tboth = array2table(tboth,'RowNames',simtext,'VariableNames',...
    {'mGlobal','sGlobal','mLC','sLC','mHCSS','sHCSS','mHCPS','sHCPS'});
writetable(Tboth,[sfile 'means_stds_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.csv'],'WriteRowNames',true);
save('means_stds_areaw_hist_aclim_zmeso200_obsglm100_strom_cope_global_biomes.mat',...
    'Tboth','tboth');

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








