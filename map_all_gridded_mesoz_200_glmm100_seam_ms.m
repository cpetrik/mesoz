% CMIP6 output 
% mesoz 200m integrations
% map all ESMs and GLMM together

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat');

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;
 
%%
figure
pcolor(gzmo_all)
shading flat
title('GF')

figure
pcolor(ozmo_all)
shading flat
title('obs')

figure
pcolor(lat_g)
shading flat
title('lat')

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
czmo_all = czmo_all * 12.01 * 1e3;
mzmo_all = mzmo_all * 12.01 * 1e3;
nzmo_all = nzmo_all * 12.01 * 1e3;
gzmo_all = gzmo_all * 12.01 * 1e3;
izmo_all = izmo_all * 12.01 * 1e3;
uzmo_all = uzmo_all * 12.01 * 1e3;
%obsglm mg C/m2 

%% Fix seam
% lon_s = lon_g;
% lat_s = lat_g;
% cz_s = czmo_all;
% mz_s = mzmo_all;
% nz_s = nzmo_all;
% gz_s = gzmo_all;
% iz_s = izmo_all;
% uz_s = uzmo_all;
% 
% lon_s(361,:) = lon_s(360,:)+1;
% lat_s(361,:) = lat_s(360,:);
% cz_s(361,:) = czmo_all(360,:)+eps;
% mz_s(361,:) = mzmo_all(360,:)+eps;
% nz_s(361,:) = nzmo_all(360,:)+eps;
% gz_s(361,:) = gzmo_all(360,:)+eps;
% iz_s(361,:) = izmo_all(360,:)+eps;
% uz_s(361,:) = uzmo_all(360,:)+eps;

[lat_s,lon_s,cz_s] = cyclic_map_seam(lat_g,lon_g,czmo_all);
[~,~,mz_s] = cyclic_map_seam(lat_g,lon_g,mzmo_all);
[~,~,nz_s] = cyclic_map_seam(lat_g,lon_g,nzmo_all);
[~,~,gz_s] = cyclic_map_seam(lat_g,lon_g,gzmo_all);
[~,~,iz_s] = cyclic_map_seam(lat_g,lon_g,izmo_all);
[~,~,uz_s] = cyclic_map_seam(lat_g,lon_g,uzmo_all);
[~,~,oz_s] = cyclic_map_seam(lat_g,lon_g,ozmo_all);

cmap = cmocean('tempo',8);

%% log trans Map
% with CMCC
% 8plot by type 
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(oz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(cz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(nz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(iz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E - 
% subplot('Position',[0.5 0.75 0.44 0.25])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(lat_s,lon_s,vp)
% cmocean('tempo')
% colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 6])
% set(gcf,'renderer','painters')
% text(0,1.75,'P','HorizontalAlignment','center')

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(mz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(gz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(uz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'UK','HorizontalAlignment','center','FontWeight','bold')

print('-dpng',[ppath 'Map_all_hist_clim_glmm100_log10mgCm-2_int200m_ms.png'])

%% colorbins
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(oz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(cz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(nz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(iz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E - 
% subplot('Position',[0.5 0.75 0.44 0.25])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(lat_s,lon_s,vp)
% colormap(cmap)
% colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 6])
% set(gcf,'renderer','painters')
% text(0,1.75,'P','HorizontalAlignment','center')

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(mz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(gz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(uz_s))
colormap(cmap)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')

print('-dpng',[ppath 'Map_all_hist_clim_glmm100_log10mgCm-2_int200m_bins.png'])

%% contours
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(oz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
contourm(lat_s,lon_s,log10(oz_s),2:4,'LineColor','k')
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(cz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(cz_s),2:4,'LineColor','k')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(nz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(nz_s),2:4,'LineColor','k')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(iz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(iz_s),2:4,'LineColor','k')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E - 
% subplot('Position',[0.5 0.75 0.44 0.25])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(lat_s,lon_s,vp)
% cmocean('tempo')
% colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 6])
% set(gcf,'renderer','painters')
% text(0,1.75,'P','HorizontalAlignment','center')

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(mz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(mz_s),2:4,'LineColor','k')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(gz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(gz_s),2:4,'LineColor','k')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,log10(uz_s))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 4])
set(gcf,'renderer','painters')
contourm(lat_s,lon_s,log10(uz_s),2:4,'LineColor','k')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')

print('-dpng',[ppath 'Map_all_hist_clim_glmm100_log10mgCm-2_int200m_contour.png'])
