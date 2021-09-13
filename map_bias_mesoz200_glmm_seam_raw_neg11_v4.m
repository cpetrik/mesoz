% CMIP6 output 
% mesoz 200m integrations
% Map ESM bias raw and scaled
% v4 = Add 1e-8 instead of eps

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_same_orientation.mat');

stex = 'add 1E-8';

%%
% %% CAN
% cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
% load([cpath 'can_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all');
% cz = zmo_all;
% cz(cz(:)<0) = 0;
% clear zmo_all 
% 
% %% CMCC
% mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';
% load([mpath 'cmcc_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all');
% mz = zmo_all;
% mz(mz(:)<0) = 0;
% clear zmo_all 
% 
% %% CNRM
% npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
% load([npath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all');
% 
% nz = zmo_all;
% nz(nz(:)<0) = 0;
% clear zmo_all 
% 
% %% UKESM
% upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
% load([upath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all');
% uz = zmo_all;
% uz(uz(:)<0) = 0;
% clear zmo_all 
% 
% %% IPSL
% ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
% load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all');
% iz = zmo_all;
% iz(iz(:)<0) = 0;
% clear zmo_all 
% 
% %% GFDL
% gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
% load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
%     'zmo_all','lat','lon');
% gz = zmo_all;
% gz(gz(:)<0) = 0;
% clear zmo_all 
% 
% %% Chl, SST, GLM zmeso, grid vars
% opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
% load([opath 'glm_obs_mesoz.mat']);
% load([opath 'glm_obs_grid.mat'])
% 
% % Reshape
% [ni,nj] = size(gz);
% glmz_mo = reshape(glmobsmesoz,ni,nj,12);
% lat_c = reshape(Lat,ni,nj);
% lon_c = reshape(Lon,ni,nj);
% 
% % annual climatology
% oz = nanmean(glmz_mo,3);
% 
% %% Convert all zoo units to mgC/m2
% %all models in molC: 12.01 g C in 1 mol C
% %1e3 mg in 1 g
% cz = cz * 12.01 * 1e3;
% mz = mz * 12.01 * 1e3;
% nz = nz * 12.01 * 1e3;
% gz = gz * 12.01 * 1e3;
% iz = iz * 12.01 * 1e3;
% uz = uz * 12.01 * 1e3;
% %obsglm mg C/m2 

%%
clatlim=[-90 90];
clonlim=[-280 80];

% [lat_g,lon_g] = meshgrid(lat,lon);
load coastlines;

%%
% figure
% pcolor(cz)
% shading flat
% title('CAN')
% 
% figure
% pcolor(mz)
% shading flat
% title('CMCC')
% 
% figure
% pcolor(nz)
% shading flat
% title('CNRM')
% 
% figure
% pcolor(iz)
% shading flat
% title('IP')
% 
% figure
% pcolor(gz)
% shading flat
% title('GF')
% 
% figure
% pcolor(uz)
% shading flat
% title('UK')
% 
% figure
% pcolor(oz)
% shading flat
% title('obs')
% 
% %% flip as needed
% close all
% cz = fliplr(cz);
% mz = fliplr(mz);
% nz = fliplr(nz);
% iz = fliplr(iz);
% uz = fliplr(uz);
% 
% %% save same orientation
% lon_o = lon_c;
% lat_o = lat_c;
% 
% save('cmip6_hist_space_means_50yr_zmeso200_same_orientation.mat',...
%     'cz','mz','nz','gz','iz','uz','oz',...
%     'lon_g','lat_g','lon_o','lat_o');

%% Fix seam
lon_s = lon_g;
lat_s = lat_g;
cz_s = cz;
mz_s = mz;
nz_s = nz;
gz_s = gz;
iz_s = iz;
uz_s = uz;

lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
cz_s(361,:) = cz(360,:)+eps;
mz_s(361,:) = mz(360,:)+eps;
nz_s(361,:) = nz(360,:)+eps;
gz_s(361,:) = gz(360,:)+eps;
iz_s(361,:) = iz(360,:)+eps;
uz_s(361,:) = uz(360,:)+eps;

oz_o = oz;
lon_o(361,:) = lon_o(360,:)+1;
lat_o(361,:) = lat_o(360,:);
oz_o(361,:) = oz(360,:)+eps;

%% Log-trans
cz_r = log(cz_s+1e-8);
mz_r = log(mz_s+1e-8);
nz_r = log(nz_s+1e-8);
gz_r = log(gz_s+1e-8);
iz_r = log(iz_s+1e-8);
uz_r = log(uz_s+1e-8);
oz_r = log(oz_o+1e-8);

%% Scale -1 to 1
Cmax = nanmax(cz_r(:));
Cmin  = nanmin(cz_r(:));
cz_n = -1 + 2.*((cz_r - Cmin) ./ (Cmax - Cmin));

Mmax = nanmax(mz_r(:));
Mmin  = nanmin(mz_r(:));
mz_n = -1 + 2.*((mz_r - Mmin) ./ (Mmax - Mmin));

Nmax = nanmax(nz_r(:));
Nmin  = nanmin(nz_r(:));
nz_n = -1 + 2.*((nz_r - Nmin) ./ (Nmax - Nmin));

Gmax = nanmax(gz_r(:));
Gmin  = nanmin(gz_r(:));
gz_n = -1 + 2.*((gz_r - Gmin) ./ (Gmax - Gmin));

Imax = nanmax(iz_r(:));
Imin  = nanmin(iz_r(:));
iz_n = -1 + 2.*((iz_r - Imin) ./ (Imax - Imin));

Umax = nanmax(uz_r(:));
Umin  = nanmin(uz_r(:));
uz_n = -1 + 2.*((uz_r - Umin) ./ (Umax - Umin));

Omax = nanmax(oz_r(:));
Omin  = nanmin(oz_r(:));
oz_n = -1 + 2.*((oz_r - Omin) ./ (Omax - Omin));

tab(1,1) = Cmin;
tab(1,2) = Cmax;
tab(2,1) = Mmin;
tab(2,2) = Mmax;
tab(3,1) = Nmin;
tab(3,2) = Nmax;
tab(4,1) = Gmin;
tab(4,2) = Gmax;
tab(5,1) = Imin;
tab(5,2) = Imax;
tab(6,1) = Umin;
tab(6,2) = Umax;
tab(7,1) = Omin;
tab(7,2) = Omax;

%% Diffs (Bias)
diff_cr = cz_r - oz_r;
diff_mr = mz_r - oz_r;
diff_nr = nz_r - oz_r;
diff_gr = gz_r - oz_r;
diff_ir = iz_r - oz_r;
diff_ur = uz_r - oz_r;

diff_cn = cz_n - oz_n;
diff_mn = mz_n - oz_n;
diff_nn = nz_n - oz_n;
diff_gn = gz_n - oz_n;
diff_in = iz_n - oz_n;
diff_un = uz_n - oz_n;

%% log trans Map
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_o,lon_o,(oz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log mgC m^-^2)')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(cz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E - 
% subplot('Position',[0.5 0.75 0.44 0.25])

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(mz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(gz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_r))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 9])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
stamp(stex)
print('-dpng',[ppath 'Map_all_hist_clim_glmm_log_scale_v4.png'])

%% -1 to 1 Map
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_o,lon_o,(oz_n))
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
surfm(lat_s,lon_s,(cz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_n))
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
surfm(lat_s,lon_s,(mz_n))
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
surfm(lat_s,lon_s,(gz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
stamp(stex)
print('-dpng',[ppath 'Map_all_hist_clim_glmm_neg11_scale_v4.png'])

%% Bias log trans Map
figure(3)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_cr)
cmocean('balance')
caxis([-2.5 2.5])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_mr)
cmocean('balance')
caxis([-2.5 2.5])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nr)
cmocean('balance')
caxis([-2.5 2.5])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gr)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'log zmeso (mgC m^-^2)')
caxis([-2.5 2.5])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ir)
cmocean('balance')
caxis([-2.5 2.5])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ur)
cmocean('balance')
caxis([-2.5 2.5])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
stamp(stex)
print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm_log_scale_v4.png'])

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
surfm(lat_s,lon_s,diff_mn)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nn)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gn)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_in)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_un)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
stamp(stex)
print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm_neg11_scale_v4.png'])

