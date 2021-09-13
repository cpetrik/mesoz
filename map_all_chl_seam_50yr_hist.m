% CMIP6 output 
% chl Hist 50yr mean

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HCchl50');

%% CMCC 
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
load([mpath 'cmcc_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HMchl50');

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HNchl50');

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HUchl50');

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HIchl50'); 

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_space_means_zmeso200_schl_sst.mat'],'HGchl50');

load([gpath 'hist/gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'lat','lon');

[lat_g,lon_g] = meshgrid(lat,lon);
clear lat lon

%% obs
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';
opath='/Volumes/MIP/Obs_Data/Chl/';
load([opath 'globcolour_soblend_mc_x1.mat'])

chl_ocx = chl*1e-3;
[lat_c,lon_c] = meshgrid(lat,lon);

mchl = nanmean(chl_ocx,3);

%% Chl, SST, GLM zmeso, grid vars
gpath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([gpath 'glm_obs_chl.mat']);
load([gpath 'glm_obs_grid.mat'])

% Reshape
[ni,nj,nt] = size(HGchl50);
oc_mo = reshape(glmobschl,ni,nj,12); 
lat_o = reshape(Lat,ni,nj);
lon_o = reshape(Lon,ni,nj);
oc = nanmean(oc_mo,3);

%% Convert all zoo units to mgC/m2
%chl in kg m-3', put in g m-3 
%chl in g m-3', put in mg m-3 (CNRM & IPSL)
%obs already in those units
HCchl = HCchl50 * 1e6;
HMchl = HMchl50 * 1e6;
HGchl = HGchl50 * 1e6;
HUchl = HUchl50 * 1e6;
HNchl = HNchl50 * 1e3;
HIchl = HIchl50 * 1e3;
mchl = mchl * 1e3;

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

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

figure %Rus R
pcolor(oc)
shading flat
title('obs glmm')

figure
pcolor(mchl)
shading flat
title('obs sat')


%% lat
close all
figure
pcolor(lat_g)
shading flat
title('latg')

figure
pcolor(lat_c)
shading flat
title('latc')

figure
pcolor(lat_o)
shading flat
title('lato')

%% flip as needed
close all
HGchl = fliplr(HGchl);
oc = fliplr(oc);
lat_g = fliplr(lat_g);

%% Fix seam
[lat_s,lon_s,Cchl] = cyclic_map_seam(lat_g,lon_g,HCchl);
[~,~,Mchl] = cyclic_map_seam(lat_g,lon_g,HMchl);
[~,~,Nchl] = cyclic_map_seam(lat_g,lon_g,HNchl);
[~,~,Gchl] = cyclic_map_seam(lat_g,lon_g,HGchl);
[~,~,Ichl] = cyclic_map_seam(lat_g,lon_g,HIchl);
[~,~,Uchl] = cyclic_map_seam(lat_g,lon_g,HUchl);
[~,~,OchlS] = cyclic_map_seam(lat_g,lon_g,mchl);
[~,~,OchlG] = cyclic_map_seam(lat_g,lon_g,oc);

%% sat chl
% 8plot by type 
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(OchlS))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'obs','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'surf chl (log_1_0 mgC m^-^3)')

%CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Cchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Nchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Ichl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Mchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Gchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Uchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_satchl_log10mgCm-3.png'])

%% glmm chl
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f3.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(OchlG))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'surf chl (log_1_0 mgC m^-^3)')

%CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Cchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Nchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Ichl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Mchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Gchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(Uchl))
cmocean('tempo')
caxis([-2 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_glmchl_log10mgCm-3.png'])

