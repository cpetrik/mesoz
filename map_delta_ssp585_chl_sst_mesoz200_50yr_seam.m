% CMIP6 output 
% mesoz 200m integrations
% Map ESM future diffs 50-yr means

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

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
HCzm = HCzm50 * 12.01 * 1e3;
HNzm = HNzm50 * 12.01 * 1e3;
HGzm = HGzm50 * 12.01 * 1e3;
HIzm = HIzm50 * 12.01 * 1e3;
HUzm = HUzm50 * 12.01 * 1e3;

FCzm = FCzm50 * 12.01 * 1e3;
FNzm = FNzm50 * 12.01 * 1e3;
FGzm = FGzm50 * 12.01 * 1e3;
FIzm = FIzm50 * 12.01 * 1e3;
FUzm = FUzm50 * 12.01 * 1e3;

%chl in kg m-3', put in g m-3 
%chl in g m-3', put in mg m-3 (CNRM & IPSL)
HCchl = HCchl50 * 1e6;
HGchl = HGchl50 * 1e6;
HUchl = HUchl50 * 1e6;
HNchl = HNchl50 * 1e3;
HIchl = HIchl50 * 1e3;

FCchl = FCchl50 * 1e6;
FGchl = FGchl50 * 1e6;
FUchl = FUchl50 * 1e6;
FNchl = FNchl50 * 1e3;
FIchl = FIchl50 * 1e3;

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
pcolor(HCsst50)
shading flat
title('CAN')

figure
pcolor(HNsst50)
shading flat
title('CNRM')

figure
pcolor(HIsst50)
shading flat
title('IP')

figure
pcolor(HGsst50)
shading flat
title('GF')

figure
pcolor(HUsst50)
shading flat
title('UK')

%% ssp sst50
close all
figure
pcolor(FCsst50)
shading flat
title('CAN')

figure
pcolor(FNsst50)
shading flat
title('CNRM')

figure
pcolor(FIsst50)
shading flat
title('IP')

figure
pcolor(FGsst50)
shading flat
title('GF')

figure
pcolor(FUsst50)
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
HIsst50 = fliplr(HIsst50);
HGsst50 = fliplr(HGsst50);
HUsst50 = fliplr(HUsst50);
FIsst50 = fliplr(FIsst50);
FGsst50 = fliplr(FGsst50);
FUsst50 = fliplr(FUsst50);
lat_g = fliplr(lat_g);

%% Diffs 
diff_ct = FCsst50 - HCsst50;
diff_nt = FNsst50 - HNsst50;
diff_gt = FGsst50 - HGsst50;
diff_it = FIsst50 - HIsst50;
diff_ut = FUsst50 - HUsst50;

%% Pdiffs
pdiff_cc = (FCchl - HCchl) ./ HCchl;
pdiff_nc = (FNchl - HNchl) ./ HNchl;
pdiff_gc = (FGchl - HGchl) ./ HGchl;
pdiff_ic = (FIchl - HIchl) ./ HIchl;
pdiff_uc = (FUchl - HUchl) ./ HUchl;

pdiff_cz = (FCzm - HCzm) ./ HCzm;
pdiff_nz = (FNzm - HNzm) ./ HNzm;
pdiff_gz = (FGzm - HGzm) ./ HGzm;
pdiff_iz = (FIzm - HIzm) ./ HIzm;
pdiff_uz = (FUzm - HUzm) ./ HUzm;

%% Fix seam
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);
[~,~,diff_nt2] = cyclic_map_seam(lat_g,lon_g,diff_nt);
[~,~,diff_gt2] = cyclic_map_seam(lat_g,lon_g,diff_gt);
[~,~,diff_it2] = cyclic_map_seam(lat_g,lon_g,diff_it);
[~,~,diff_ut2] = cyclic_map_seam(lat_g,lon_g,diff_ut);

[~,~,pdiff_cc2] = cyclic_map_seam(lat_g,lon_g,pdiff_cc);
[~,~,pdiff_nc2] = cyclic_map_seam(lat_g,lon_g,pdiff_nc);
[~,~,pdiff_gc2] = cyclic_map_seam(lat_g,lon_g,pdiff_gc);
[~,~,pdiff_ic2] = cyclic_map_seam(lat_g,lon_g,pdiff_ic);
[~,~,pdiff_uc2] = cyclic_map_seam(lat_g,lon_g,pdiff_uc);

[~,~,pdiff_cz2] = cyclic_map_seam(lat_g,lon_g,pdiff_cz);
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
caxis([-8 8])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nt2)
cmocean('balance')
caxis([-8 8])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gt2)
cmocean('balance')
caxis([-8 8])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_it2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'SST (^oC)')
caxis([-8 8])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ut2)
cmocean('balance')
caxis([-8 8])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_SST_diff_Hist_SSP585_50yr.png'])

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
surfm(lat_s,lon_s,pdiff_nc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_ic2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'% \Delta chl (mg m^-^3)')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uc2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_chl_pdiff_Hist_SSP585_50yr.png'])

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
surfm(lat_s,lon_s,pdiff_nz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'% \Delta mesoz (mgC m^-^2)')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_zmeso_pdiff_Hist_SSP585_50yr.png'])

