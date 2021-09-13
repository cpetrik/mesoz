% Calculate different skill metrics for each ESM
% log transform biomass
% uses SH shift

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load('skill_hist_model_obsglm_climatols.mat')
%load('skill_model_obsglm_climatols.mat')

%% Standardization
% 1st take log
Acomb = log(comb(:,3:9) + eps);
Dcomb = log(dcomb(:,3:9) + eps);
Jcomb = log(jcomb(:,3:9) + eps);
Mcomb = log(mcomb(:,3:9) + eps);
Scomb = log(scomb(:,3:9) + eps);

% then stdize to -1 to 1
%sc = -1 + 2.*(data - min(data))./(max(data) - min(data));
Amax = nanmax(Acomb);
Amin  = nanmin(Acomb);
Acomb = -1 + 2.*((Acomb - Amin) ./ (Amax - Amin));

%%
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'lat','lon')
[lat_g,lon_g] = meshgrid(lat,lon);

%%
cz = Acomb(:,1);
mz = Acomb(:,2);
nz = Acomb(:,3);
gz = Acomb(:,4);
iz = Acomb(:,5);
uz = Acomb(:,6);
oz = Acomb(:,7);

%% Reshape
[ni,nj] = size(lon_g);
cz = reshape(cz,ni,nj);
mz = reshape(mz,ni,nj);
nz = reshape(nz,ni,nj);
gz = reshape(gz,ni,nj);
iz = reshape(iz,ni,nj);
uz = reshape(uz,ni,nj);
oz = reshape(oz,ni,nj);

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

%%
figure
pcolor(cz)
shading flat
title('CAN')

figure
pcolor(mz)
shading flat
title('CMCC')

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

%% Fix seam
close all

lon_s = lon_g;
lat_s = lat_g;
cz_s = cz;
mz_s = mz;
nz_s = nz;
gz_s = gz;
iz_s = iz;
uz_s = uz;
oz_s = oz;

lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
cz_s(361,:) = cz(360,:)+eps;
mz_s(361,:) = mz(360,:)+eps;
nz_s(361,:) = nz(360,:)+eps;
gz_s(361,:) = gz(360,:)+eps;
iz_s(361,:) = iz(360,:)+eps;
uz_s(361,:) = uz(360,:)+eps;
oz_s(361,:) = oz(360,:)+eps;

%% Diffs (Bias)
diff_cn = cz_s - oz_s;
diff_mn = mz_s - oz_s;
diff_nn = nz_s - oz_s;
diff_gn = gz_s - oz_s;
diff_in = iz_s - oz_s;
diff_un = uz_s - oz_s;

%% -1 to 1 Map
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(oz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(cz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E - 
% subplot('Position',[0.5 0.75 0.44 0.25])

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(mz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(gz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
print('-dpng',[ppath 'Map_all_hist_clim_glmm_neg11_scale_v2.png'])

%% Bias -1 to 1 Map
figure(4)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_cn)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nn)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gn)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_in)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_un)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_un)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm_neg11_scale_v2.png'])

