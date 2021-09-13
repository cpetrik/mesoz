% CMIP6 output 
% 200m integrations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','units');

cz = zmo_all;
cz(cz(:)<0) = 0;
cunits = units;

clear zmo_all units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','units');

nz = zmo_all;
nz(nz(:)<0) = 0;
nunits = units;

clear zmo_all units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','units');

uz = zmo_all;
uz(uz(:)<0) = 0;
%uunits = units;

clear zmo_all units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','units');

iz = zmo_all;
iz(iz(:)<0) = 0;
%iunits = units;

clear zmo_all units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','lat','lon');

gz = zmo_all;
gz(gz(:)<0) = 0;

clear zmo_all 

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_mesoz.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj] = size(gz);
glmz_mo = reshape(glmobsmesoz,ni,nj,12);
lat_c = reshape(Lat,ni,nj);
lon_c = reshape(Lon,ni,nj);

% annual climatology
oz = nanmean(glmz_mo,3);

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
cz = cz * 12.01 * 1e3;
nz = nz * 12.01 * 1e3;
gz = gz * 12.01 * 1e3;
iz = iz * 12.01 * 1e3;
uz = uz * 12.01 * 1e3;
%obsglm mg C/m2 

%%
clatlim=[-90 90];
clonlim=[-280 80];

[lat_g,lon_g] = meshgrid(lat,lon);
load coastlines;

%% Maps
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,fliplr(cz))
cmocean('tempo')
caxis([0 5e3])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,fliplr(nz))
cmocean('tempo')
caxis([0 5e3])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(gz))
cmocean('tempo')
caxis([0 5e3])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,fliplr(iz))
cmocean('tempo')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (mgC m^-^2)')
caxis([0 5e3])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,fliplr(uz))
cmocean('tempo')
caxis([0 5e3])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_c,lon_c,(oz))
cmocean('tempo')
caxis([0 5e3])
text(0.2,1.65,'obs-GLMM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_clim_glmm_mgCm-2_int200m.png'])


%% log trans
figure(2)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(fliplr(cz)))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(fliplr(nz)))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(gz))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(fliplr(iz)))
cmocean('tempo')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')
caxis([2 4])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(fliplr(uz)))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_c,lon_c,log10(oz))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'obs-GLMM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_clim_glmm_log10mgCm-2_int200m.png'])

%% Test getting rid of seam
lon_s = lon_g;
lat_s = lat_g;
gz_s = gz;
lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
gz_s(361,:) = gz(360,:)+eps;

%%
figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(gz_s))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
