% CMIP6 output 
% mesoz 200m integrations
% Map ESM future diffs 2090s vs. 1990s

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
load([mpath 'cmcc_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_space_means_zmeso200_schl_sst.mat']); 

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

load([gpath 'hist/gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'lat','lon');

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
HCzm = HCzm90 * 12.01 * 1e3;
HMzm = HMzm90 * 12.01 * 1e3;
HNzm = HNzm90 * 12.01 * 1e3;
HGzm = HGzm90 * 12.01 * 1e3;
HIzm = HIzm90 * 12.01 * 1e3;
HUzm = HUzm90 * 12.01 * 1e3;

FCzm = FCzm90 * 12.01 * 1e3;
FMzm = FMzm90 * 12.01 * 1e3;
FNzm = FNzm90 * 12.01 * 1e3;
FGzm = FGzm90 * 12.01 * 1e3;
FIzm = FIzm90 * 12.01 * 1e3;
FUzm = FUzm90 * 12.01 * 1e3;

%chl in kg m-3', put in g m-3 
%chl in g m-3', put in mg m-3 (CNRM & IPSL)
HCchl = HCchl90 * 1e6;
HMchl = HMchl90 * 1e6;
HGchl = HGchl90 * 1e6;
HUchl = HUchl90 * 1e6;
HNchl = HNchl90 * 1e3;
HIchl = HIchl90 * 1e3;

FCchl = FCchl90 * 1e6;
FMchl = FMchl90 * 1e6;
FGchl = FGchl90 * 1e6;
FUchl = FUchl90 * 1e6;
FNchl = FNchl90 * 1e3;
FIchl = FIchl90 * 1e3;

%%
clatlim=[-90 90];
clonlim=[-280 80];

[lat_g,lon_g] = meshgrid(lat,lon);
load coastlines;

%% hist z
figure
pcolor(HCzm)
shading flat
title('CAN')

figure
pcolor(HMzm)
shading flat
title('CMCC')

figure
pcolor(HNzm)
shading flat
title('CNRM')

figure
pcolor(HIzm)
shading flat
title('IP')

figure
pcolor(HGzm)
shading flat
title('GF')

figure
pcolor(HUzm)
shading flat
title('UK')

%% ssp z
close all
figure
pcolor(FCzm)
shading flat
title('CAN')

figure
pcolor(FMzm)
shading flat
title('CMCC')

figure
pcolor(FNzm)
shading flat
title('CNRM')

figure
pcolor(FIzm)
shading flat
title('IP')

figure
pcolor(FGzm)
shading flat
title('GF')

figure
pcolor(FUzm)
shading flat
title('UK')

%% hist chl
close all
figure
pcolor(HCchl)
shading flat
title('CAN')

figure
pcolor(HMchl)
shading flat
title('CMCC')

figure
pcolor(HNchl)
shading flat
title('CNRM')

figure
pcolor(HIchl)
shading flat
title('IP')

figure
pcolor(HGchl)
shading flat
title('GF')

figure
pcolor(HUchl)
shading flat
title('UK')

%% ssp chl
close all
figure
pcolor(FCchl)
shading flat
title('CAN')

figure
pcolor(FMchl)
shading flat
title('CMCC')

figure
pcolor(FNchl)
shading flat
title('CNRM')

figure
pcolor(FIchl)
shading flat
title('IP')

figure
pcolor(FGchl)
shading flat
title('GF')

figure
pcolor(FUchl)
shading flat
title('UK')

%% hist sst
close all
figure
pcolor(HCsst90)
shading flat
title('CAN')

figure
pcolor(HMsst90)
shading flat
title('CMCC')

figure
pcolor(HNsst90)
shading flat
title('CNRM')

figure
pcolor(HIsst90)
shading flat
title('IP')

figure
pcolor(HGsst90)
shading flat
title('GF')

figure
pcolor(HUsst90)
shading flat
title('UK')

%% ssp sst90
close all
figure
pcolor(FCsst90)
shading flat
title('CAN')

figure
pcolor(FMsst90)
shading flat
title('CMCC')

figure
pcolor(FNsst90)
shading flat
title('CNRM')

figure
pcolor(FIsst90)
shading flat
title('IP')

figure
pcolor(FGsst90)
shading flat
title('GF')

figure
pcolor(FUsst90)
shading flat
title('UK')

%% lat
close all
figure
pcolor(lat_g)
shading flat
title('lat')

%% flip as needed
close all
HGzm = fliplr(HGzm);
FGzm = fliplr(FGzm);
FUzm = fliplr(FUzm);
HGchl = fliplr(HGchl);
FGchl = fliplr(FGchl);
HIsst90 = fliplr(HIsst90);
HGsst90 = fliplr(HGsst90);
HUsst90 = fliplr(HUsst90);
FIsst90 = fliplr(FIsst90);
FGsst90 = fliplr(FGsst90);
FUsst90 = fliplr(FUsst90);
lat_g = fliplr(lat_g);

%% save same orientation
save('cmip6_hist_ssp585_space_means_90s_zmeso200_schl_sst_same_orientation.mat',...
    'HCzm','HMzm','HNzm','HGzm','HIzm','HUzm',...
    'FCzm','FMzm','FNzm','FGzm','FIzm','FUzm',...
    'HCchl','HMchl','HNchl','HGchl','HIchl','HUchl',...
    'FCchl','FMchl','FNchl','FGchl','FIchl','FUchl',...
    'HCsst90','HMsst90','HNsst90','HGsst90','HIsst90','HUsst90',...
    'FCsst90','FMsst90','FNsst90','FGsst90','FIsst90','FUsst90',...
    'lat_g','lon_g');

%% Diffs 
diff_ct = FCsst90 - HCsst90;
diff_mt = FMsst90 - HMsst90;
diff_nt = FNsst90 - HNsst90;
diff_gt = FGsst90 - HGsst90;
diff_it = FIsst90 - HIsst90;
diff_ut = FUsst90 - HUsst90;

%% Pdiffs
pdiff_cc = (FCchl - HCchl) ./ HCchl;
pdiff_mc = (FMchl - HMchl) ./ HMchl;
pdiff_nc = (FNchl - HNchl) ./ HNchl;
pdiff_gc = (FGchl - HGchl) ./ HGchl;
pdiff_ic = (FIchl - HIchl) ./ HIchl;
pdiff_uc = (FUchl - HUchl) ./ HUchl;

pdiff_cz = (FCzm - HCzm) ./ HCzm;
pdiff_mz = (FMzm - HMzm) ./ HMzm;
pdiff_nz = (FNzm - HNzm) ./ HNzm;
pdiff_gz = (FGzm - HGzm) ./ HGzm;
pdiff_iz = (FIzm - HIzm) ./ HIzm;
pdiff_uz = (FUzm - HUzm) ./ HUzm;

%% Fix seam
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);
[~,~,diff_mt2] = cyclic_map_seam(lat_g,lon_g,diff_mt);
[~,~,diff_nt2] = cyclic_map_seam(lat_g,lon_g,diff_nt);
[~,~,diff_gt2] = cyclic_map_seam(lat_g,lon_g,diff_gt);
[~,~,diff_it2] = cyclic_map_seam(lat_g,lon_g,diff_it);
[~,~,diff_ut2] = cyclic_map_seam(lat_g,lon_g,diff_ut);

[~,~,pdiff_cc2] = cyclic_map_seam(lat_g,lon_g,pdiff_cc);
[~,~,pdiff_mc2] = cyclic_map_seam(lat_g,lon_g,pdiff_mc);
[~,~,pdiff_nc2] = cyclic_map_seam(lat_g,lon_g,pdiff_nc);
[~,~,pdiff_gc2] = cyclic_map_seam(lat_g,lon_g,pdiff_gc);
[~,~,pdiff_ic2] = cyclic_map_seam(lat_g,lon_g,pdiff_ic);
[~,~,pdiff_uc2] = cyclic_map_seam(lat_g,lon_g,pdiff_uc);

[~,~,pdiff_cz2] = cyclic_map_seam(lat_g,lon_g,pdiff_cz);
[~,~,pdiff_mz2] = cyclic_map_seam(lat_g,lon_g,pdiff_mz);
[~,~,pdiff_nz2] = cyclic_map_seam(lat_g,lon_g,pdiff_nz);
[~,~,pdiff_gz2] = cyclic_map_seam(lat_g,lon_g,pdiff_gz);
[~,~,pdiff_iz2] = cyclic_map_seam(lat_g,lon_g,pdiff_iz);
[~,~,pdiff_uz2] = cyclic_map_seam(lat_g,lon_g,pdiff_uz);

%% SST diff Map
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ct2)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_mt2)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nt2)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gt2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'SST (^oC)')
caxis([-10 10])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_it2)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ut2)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_SST_diff_Hist_SSP585_90s.png'])

%% Chl pdiff Map
figure(2)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gc2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'% \Delta chl (mg m^-^3)')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_ic2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_chl_pdiff_Hist_SSP585_90s.png'])

%% Zmeso pdiff Map
figure(3)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'% \Delta mesoz (mgC m^-^2)')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_zmeso_pdiff_Hist_SSP585_90s.png'])

%% put chl and zoo next to each other
f5 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Chl cmcc
subplot('Position',[0.025 0.8 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mc2)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'% \Delta chl (mg m^-^3)','HorizontalAlignment','center','FontWeight','bold')
text(-1.95,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Chl cnrm
subplot('Position',[0.025 0.6 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nc2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Chl gfdl
subplot('Position',[0.025 0.4 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gc2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Chl IPSL
subplot('Position',[0.025 0.2 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_ic2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Chl uk
subplot('Position',[0.025 0.0 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uc2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - Zoo CMCC
subplot('Position',[0.43 0.8 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mz2)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'% \Delta mesoz (mgC m^-^2)','HorizontalAlignment','center','FontWeight','bold')
text(-1.95,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%7 - Zoo CNRM
subplot('Position',[0.43 0.6 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nz2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Zoo GFDL
subplot('Position',[0.43 0.4 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.85 0.3 0.025 0.4]);

%9 - Zoo IPSL
subplot('Position',[0.43 0.2 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%10 - Zoo uk
subplot('Position',[0.43 0.0 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_chl_zmeso_pdiff_Hist_SSP585_90s_cmcc.png'])

%% Just CAN for Supp
figure(8)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cc2)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'% \Delta chl (mg m^-^3)','HorizontalAlignment','center','FontWeight','bold')
%text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cz2)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'% \Delta mesoz (mgC m^-^2)','HorizontalAlignment','center','FontWeight','bold')
%text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.8 0.68 0.025 0.3]);
print('-dpng',[ppath 'Map_chl_zmeso_pdiff_Hist_SSP585_90s_CANonly.png'])


%% some sort of trophic amp calc
f6 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Chl can
subplot('Position',[0.025 0.8 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cz2./pdiff_cc2)
cmocean('balance')
caxis([-4 4])
text(0,1.75,'Ratio mesoz:chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Chl cnrm
subplot('Position',[0.025 0.6 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nz2./pdiff_nc2)
cmocean('balance')
caxis([-4 4])
text(-1.85,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Chl gfdl
subplot('Position',[0.025 0.4 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2./pdiff_gc2)
cmocean('balance')
caxis([-4 4])
text(-1.85,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cbr = colorbar('Position',[0.85 0.7 0.025 0.25]);
xlabel(cbr,'ratio')

%4 - Chl IPSL
subplot('Position',[0.025 0.2 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2./pdiff_ic2)
cmocean('balance')
caxis([-4 4])
text(-1.85,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Chl uk
subplot('Position',[0.025 0.0 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2./pdiff_uc2)
cmocean('balance')
caxis([-4 4])
text(-1.85,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - Zoo CAN
subplot('Position',[0.43 0.8 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cz2-pdiff_cc2)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'mesoz - chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%7 - Zoo CNRM
subplot('Position',[0.43 0.6 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nz2-pdiff_nc2)
cmocean('balance')
caxis([-1 1])
text(-1.85,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Zoo GFDL
subplot('Position',[0.43 0.4 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2-pdiff_gc2)
cmocean('balance')
caxis([-1 1])
text(-1.85,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cbd = colorbar('Position',[0.85 0.1 0.025 0.25]);
xlabel(cbd,'difference')

%9 - Zoo IPSL
subplot('Position',[0.43 0.2 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2-pdiff_ic2)
cmocean('balance')
caxis([-1 1])
text(-1.85,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%10 - Zoo uk
subplot('Position',[0.43 0.0 0.4 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2-pdiff_uc2)
cmocean('balance')
caxis([-1 1])
text(-1.85,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_chl_zmeso_troph_amp_Hist_SSP585_90s.png'])
