% Change in SSP8.5 vs Hist by biome
% 2051-2100 vs. 1966-2015
% try both Hist and RCP biomes

clear all
close all

%% load Hist and SSP zmeso 
load('cmip6_hist_ssp585_space_means_50yr_zmeso200_schl_sst_same_orientation.mat',...
    'HCzm','HMzm','HNzm','HGzm','HIzm','HUzm',...
    'FCzm','FMzm','FNzm','FGzm','FIzm','FUzm',...
    'lat_g','lon_g')

%% load obs-glmm climatol
load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat',...
    'ozmo_all');

%% load Stromberg
load('Hist_mesoz_clim_biomes_samegrid.mat','sz_mo','sbiomes');
szmo_all = nanmean(sz_mo,3);
szmo_all = fliplr(szmo_all);
sbiomes = fliplr(sbiomes);

%% Hist biomes
load('Hist_mesoz_clim_biomes_samegrid.mat','hcbiomes','hmbiomes','hnbiomes',...
    'hgbiomes','hibiomes','hubiomes','biomes')

%% load SSP biomes
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/ssp585/';

load([fpath 'CanESM_ssp585_biomes.mat'],'scbiomes');
load([mpath 'CMCC_ssp585_biomes.mat'],'smbiomes');
load([fpath 'CNRM_ssp585_biomes.mat'],'snbiomes');
load([fpath 'GFDL_ssp585_biomes.mat'],'sgbiomes');
load([fpath 'IPSL_ssp585_biomes.mat'],'sibiomes');
load([fpath 'UKESM_ssp585_biomes.mat'],'subiomes','lat','lon');

[lat_b,lon_b] = meshgrid(lat,lon);
clear lat lon

%% check orientations
% hist biomes
close all
figure(1)
pcolor(hcbiomes); shading flat;
title('CAN')

figure(7)
pcolor(hmbiomes); shading flat;
title('CMCC')

figure(2)
pcolor(hnbiomes); shading flat;
title('CNRM')

figure(3)
pcolor(hibiomes); shading flat;
title('IPSL')

figure(4)
pcolor(hgbiomes); shading flat;
title('GFDL')

figure(5)
pcolor(hubiomes); shading flat;
title('UK')

figure(6)
pcolor(biomes); shading flat;
title('obs biome')

figure(7)
pcolor(ozmo_all); shading flat;
title('obsglmm')

%% ssp biomes
close all
figure(1)
pcolor(scbiomes); shading flat;
title('CAN')

figure(7)
pcolor(smbiomes); shading flat;
title('CMCC')

figure(2)
pcolor(snbiomes); shading flat;
title('CNRM')

figure(3)
pcolor(sibiomes); shading flat;
title('IPSL')

figure(4)
pcolor(sgbiomes); shading flat;
title('GFDL')

figure(5)
pcolor(subiomes); shading flat;
title('UK')

figure(6)
pcolor(lat_b); shading flat;
title('latb')

figure(7)
pcolor(lat_g); shading flat;
title('latg')

%%
ozmo_all = fliplr(ozmo_all);
biomes = fliplr(biomes);
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

%% Save SSP biomes same orient
close all
save('Hist_SSP585_biomes_samegrid.mat','hcbiomes','hmbiomes','hnbiomes',...
    'hgbiomes','hibiomes','hubiomes','biomes','scbiomes','smbiomes','snbiomes',...
    'sgbiomes','sibiomes','subiomes','sbiomes','lat_g','lon_g')

%% Diffs and Pdiffs
% Area-weighted means (only matters for diff b/c cancels out of pdiff)
apath='/Volumes/MIP/Fish-MIP/CMIP6/';
load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

diff_cz = cell_area.*(FCzm - HCzm);
diff_mz = cell_area.*(FMzm - HMzm);
diff_nz = cell_area.*(FNzm - HNzm);
diff_gz = cell_area.*(FGzm - HGzm);
diff_iz = cell_area.*(FIzm - HIzm);
diff_uz = cell_area.*(FUzm - HUzm);

pdiff_cz = (FCzm - HCzm) ./ HCzm;
pdiff_mz = (FMzm - HMzm) ./ HMzm;
pdiff_nz = (FNzm - HNzm) ./ HNzm;
pdiff_gz = (FGzm - HGzm) ./ HGzm;
pdiff_iz = (FIzm - HIzm) ./ HIzm;
pdiff_uz = (FUzm - HUzm) ./ HUzm;

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

%% Isolate Diffs
Czlc = diff_cz(cid);
Mzlc = diff_mz(mid);
Nzlc = diff_nz(nid);
Gzlc = diff_gz(gid);
Izlc = diff_iz(iid);
Uzlc = diff_uz(uid);

Calc = (cell_area(cid));
Malc = (cell_area(mid));
Nalc = (cell_area(nid));
Galc = (cell_area(gid));
Ialc = (cell_area(iid));
Ualc = (cell_area(uid));

Czss = (diff_cz(css) );
Mzss = (diff_mz(mss) );
Nzss = (diff_nz(nss) );
Gzss = (diff_gz(gss) );
Izss = (diff_iz(iss) );
Uzss = (diff_uz(uss) );

Cass = (cell_area(css));
Mass = (cell_area(mss));
Nass = (cell_area(nss));
Gass = (cell_area(gss));
Iass = (cell_area(iss));
Uass = (cell_area(uss));

Czps = (diff_cz(cps) );
Mzps = (diff_mz(mps) );
Nzps = (diff_nz(nps) );
Gzps = (diff_gz(gps) );
Izps = (diff_iz(ips) );
Uzps = (diff_uz(ups) );

Caps = (cell_area(cps));
Maps = (cell_area(mps));
Naps = (cell_area(nps));
Gaps = (cell_area(gps));
Iaps = (cell_area(ips));
Uaps = (cell_area(ups));

%% Take area-weighted means of diffs by biome
zmeans_diff(1,1) = nansum(diff_cz(:)) ./ nansum(cell_area(:));
zmeans_diff(2,1) = nansum(diff_mz(:)) ./ nansum(cell_area(:));
zmeans_diff(3,1) = nansum(diff_nz(:)) ./ nansum(cell_area(:));
zmeans_diff(4,1) = nansum(diff_gz(:)) ./ nansum(cell_area(:));
zmeans_diff(5,1) = nansum(diff_iz(:)) ./ nansum(cell_area(:));
zmeans_diff(6,1) = nansum(diff_uz(:)) ./ nansum(cell_area(:));

% LC
zmeans_diff(1,2) = nansum(Czlc(:)) ./ nansum(Calc(:));
zmeans_diff(2,2) = nansum(Mzlc(:)) ./ nansum(Malc(:));
zmeans_diff(3,2) = nansum(Nzlc(:)) ./ nansum(Nalc(:));
zmeans_diff(4,2) = nansum(Gzlc(:)) ./ nansum(Galc(:));
zmeans_diff(5,2) = nansum(Izlc(:)) ./ nansum(Ialc(:));
zmeans_diff(6,2) = nansum(Uzlc(:)) ./ nansum(Ualc(:));

% SS
zmeans_diff(1,3) = nansum(Czss(:)) ./ nansum(Cass(:));
zmeans_diff(2,3) = nansum(Mzss(:)) ./ nansum(Mass(:));
zmeans_diff(3,3) = nansum(Nzss(:)) ./ nansum(Nass(:));
zmeans_diff(4,3) = nansum(Gzss(:)) ./ nansum(Gass(:));
zmeans_diff(5,3) = nansum(Izss(:)) ./ nansum(Iass(:));
zmeans_diff(6,3) = nansum(Uzss(:)) ./ nansum(Uass(:));

% PS
zmeans_diff(1,4) = nansum(Czps(:)) ./ nansum(Caps(:));
zmeans_diff(2,4) = nansum(Mzps(:)) ./ nansum(Maps(:));
zmeans_diff(3,4) = nansum(Nzps(:)) ./ nansum(Naps(:));
zmeans_diff(4,4) = nansum(Gzps(:)) ./ nansum(Gaps(:));
zmeans_diff(5,4) = nansum(Izps(:)) ./ nansum(Iaps(:));
zmeans_diff(6,4) = nansum(Uzps(:)) ./ nansum(Uaps(:));

%%
simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK'};
sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

Tmeanz = array2table(zmeans_diff,'RowNames',simtext,'VariableNames',...
    {'Global','LC','HCSS','HCPS'});
writetable(Tmeanz,[sfile 'delta_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes.csv'],'WriteRowNames',true);
save('delta_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes.mat',...
    'Tmeanz','zmeans_diff');

%% Isolate Pdiffs
Czlc = pdiff_cz(cid);
Mzlc = pdiff_mz(mid);
Nzlc = pdiff_nz(nid);
Gzlc = pdiff_gz(gid);
Izlc = pdiff_iz(iid);
Uzlc = pdiff_uz(uid);

Czss = (pdiff_cz(css) );
Mzss = (pdiff_mz(mss) );
Nzss = (pdiff_nz(nss) );
Gzss = (pdiff_gz(gss) );
Izss = (pdiff_iz(iss) );
Uzss = (pdiff_uz(uss) );

Czps = (pdiff_cz(cps) );
Mzps = (pdiff_mz(mps) );
Nzps = (pdiff_nz(nps) );
Gzps = (pdiff_gz(gps) );
Izps = (pdiff_iz(ips) );
Uzps = (pdiff_uz(ups) );

%% Take means by biome
zmeans_pdiff(1,1) = nanmean(pdiff_cz(:));
zmeans_pdiff(2,1) = nanmean(pdiff_mz(:));
zmeans_pdiff(3,1) = nanmean(pdiff_nz(:));
zmeans_pdiff(4,1) = nanmean(pdiff_gz(:));
zmeans_pdiff(5,1) = nanmean(pdiff_iz(:));
zmeans_pdiff(6,1) = nanmean(pdiff_uz(:));

% LC
zmeans_pdiff(1,2) = nanmean(Czlc(:));
zmeans_pdiff(2,2) = nanmean(Mzlc(:));
zmeans_pdiff(3,2) = nanmean(Nzlc(:));
zmeans_pdiff(4,2) = nanmean(Gzlc(:));
zmeans_pdiff(5,2) = nanmean(Izlc(:));
zmeans_pdiff(6,2) = nanmean(Uzlc(:));

% SS
zmeans_pdiff(1,3) = nanmean(Czss(:));
zmeans_pdiff(2,3) = nanmean(Mzss(:));
zmeans_pdiff(3,3) = nanmean(Nzss(:));
zmeans_pdiff(4,3) = nanmean(Gzss(:));
zmeans_pdiff(5,3) = nanmean(Izss(:));
zmeans_pdiff(6,3) = nanmean(Uzss(:));

% PS
zmeans_pdiff(1,4) = nanmean(Czps(:)) ;
zmeans_pdiff(2,4) = nanmean(Mzps(:)) ;
zmeans_pdiff(3,4) = nanmean(Nzps(:)) ;
zmeans_pdiff(4,4) = nanmean(Gzps(:)) ;
zmeans_pdiff(5,4) = nanmean(Izps(:)) ;
zmeans_pdiff(6,4) = nanmean(Uzps(:)) ;




