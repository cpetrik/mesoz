% CMIP6 output 
% mesoz 200m integrations
% map all ESMs and GLMM together

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

clear zmo_all 

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','units');

iz = zmo_all;
iz(iz(:)<0) = 0;
%iunits = units;

clear zmo_all 

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

%%
figure
pcolor(cz)
shading flat
title('CAN')

figure
pcolor(nz)
shading flat
title('CNRM')

figure
pcolor(iz)
shading flat
title('IP')

figure
pcolor(gz)
shading flat
title('GF')

figure
pcolor(uz)
shading flat
title('UK')

figure
pcolor(oz)
shading flat
title('obs')

%% flip as needed
close all
cz = fliplr(cz);
nz = fliplr(nz);
iz = fliplr(iz);
uz = fliplr(uz);

%% Fix seam
lon_s = lon_g;
lat_s = lat_g;
cz_s = cz;
nz_s = nz;
gz_s = gz;
iz_s = iz;
uz_s = uz;

lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
cz_s(361,:) = cz(360,:)+eps;
nz_s(361,:) = nz(360,:)+eps;
gz_s(361,:) = gz(360,:)+eps;
iz_s(361,:) = iz(360,:)+eps;
uz_s(361,:) = uz(360,:)+eps;

lon_o = lon_c;
lat_o = lat_c;
oz_o = oz;
lon_o(361,:) = lon_c(360,:)+1;
lat_o(361,:) = lat_c(360,:);
oz_o(361,:) = oz(360,:)+eps;


%% log trans Map
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(cz_s))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(nz_s))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(gz_s))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(iz_s))
cmocean('tempo')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')
caxis([2 4])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,log10(uz_s))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_o,lon_o,log10(oz_o))
cmocean('tempo')
caxis([2 4])
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_clim_glmm_log10mgCm-2_int200m.png'])

