% CMIP6 output 
% mesoz 200m integrations
% Map ESM bias norm

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_same_orientation.mat');

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
cz_r = log(cz_s+eps);
mz_r = log(mz_s+eps);
nz_r = log(nz_s+eps);
gz_r = log(gz_s+eps);
iz_r = log(iz_s+eps);
uz_r = log(uz_s+eps);
oz_r = log(oz_o+eps);

%% Scale bt normalizing
% Amean = nanmean(Acomb);
% Astd  = nanstd(Acomb);
% Acomb = (Acomb - Amean) ./ Astd;

Cmean = nanmean(cz_r(:));
Cstd  = nanstd(cz_r(:));
cz_n = (cz_r - Cmean) ./ (Cstd);

Mmean = nanmean(mz_r(:));
Mstd  = nanstd(mz_r(:));
mz_n = (mz_r - Mstd) ./ (Mstd);

Nmean = nanmean(nz_r(:));
Nstd  = nanstd(nz_r(:));
nz_n = (nz_r - Nstd) ./ (Nstd);

Gmean = nanmean(gz_r(:));
Gstd  = nanstd(gz_r(:));
gz_n = (gz_r - Gstd) ./ (Gstd);

Imean = nanmean(iz_r(:));
Istd  = nanstd(iz_r(:));
iz_n = (iz_r - Istd) ./ (Istd);

Umean = nanmean(uz_r(:));
Ustd  = nanstd(uz_r(:));
uz_n = (uz_r - Ustd) ./ (Ustd);

Omean = nanmean(oz_r(:));
Ostd  = nanstd(oz_r(:));
oz_n = (oz_r - Ostd) ./ (Ostd);

tab(1,1) = Cstd;
tab(1,2) = Cmean;
tab(2,1) = Mstd;
tab(2,2) = Mmean;
tab(3,1) = Nstd;
tab(3,2) = Nmean;
tab(4,1) = Gstd;
tab(4,2) = Gmean;
tab(5,1) = Istd;
tab(5,2) = Imean;
tab(6,1) = Ustd;
tab(6,2) = Umean;
tab(7,1) = Ostd;
tab(7,2) = Omean;

%% Diffs (Bias)
diff_cn = cz_n - oz_n;
diff_mn = mz_n - oz_n;
diff_nn = nz_n - oz_n;
diff_gn = gz_n - oz_n;
diff_in = iz_n - oz_n;
diff_un = uz_n - oz_n;

%% Norm Map
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_o,lon_o,(oz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 15])
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
caxis([1 15])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 15])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 15])
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
caxis([1 15])
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
caxis([1 15])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_n))
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 15])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
print('-dpng',[ppath 'Map_all_hist_clim_glmm_norm_scale.png'])

%% Bias norm Map
figure(4)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_cn)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_mn)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nn)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gn)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso')
caxis([-10 10])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_in)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_un)
cmocean('balance')
caxis([-10 10])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm_norm_scale.png'])

