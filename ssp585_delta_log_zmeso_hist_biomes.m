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

%%  biomes
load('Hist_SSP585_biomes_samegrid.mat','hcbiomes','hmbiomes','hnbiomes',...
    'hgbiomes','hibiomes','hubiomes')

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
pcolor(HMzm); shading flat;
title('Mz')

figure(7)
pcolor(lat_g); shading flat;
title('latg')

%% Diffs and Pdiffs
% Area-weighted means (only matters for diff b/c cancels out of pdiff)
apath='/Volumes/MIP/Fish-MIP/CMIP6/';
load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

diff_cz = cell_area.*(log10(FCzm) - log10(HCzm));
diff_mz = cell_area.*(log10(FMzm) - log10(HMzm));
diff_nz = cell_area.*(log10(FNzm) - log10(HNzm));
diff_gz = cell_area.*(log10(FGzm) - log10(HGzm));
diff_iz = cell_area.*(log10(FIzm) - log10(HIzm));
diff_uz = cell_area.*(log10(FUzm) - log10(HUzm));

% pdiff_cz = (FCzm - HCzm) ./ HCzm;
% pdiff_mz = (FMzm - HMzm) ./ HMzm;
% pdiff_nz = (FNzm - HNzm) ./ HNzm;
% pdiff_gz = (FGzm - HGzm) ./ HGzm;
% pdiff_iz = (FIzm - HIzm) ./ HIzm;
% pdiff_uz = (FUzm - HUzm) ./ HUzm;

%% ID each biome
% Exclude LC > 45N/S
pol = find(lat_g <= 45 & lat_g >= -45);

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
writetable(Tmeanz,[sfile 'delta_log10_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes.csv'],'WriteRowNames',true);
save('delta_log10_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes.mat',...
    'Tmeanz','zmeans_diff');

