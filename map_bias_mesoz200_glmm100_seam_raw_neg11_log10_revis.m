% CMIP6 output
% mesoz 200m integrations
% Map ESM bias raw and scaled
% log-10 trans

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat');

stex = 'log10';

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

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

%% fix seam
[lat_s,lon_s,cz_s] = cyclic_map_seam(lat_g,lon_g,czmo_all);
[~,~,mz_s] = cyclic_map_seam(lat_g,lon_g,mzmo_all);
[~,~,nz_s] = cyclic_map_seam(lat_g,lon_g,nzmo_all);
[~,~,gz_s] = cyclic_map_seam(lat_g,lon_g,gzmo_all);
[~,~,iz_s] = cyclic_map_seam(lat_g,lon_g,izmo_all);
[~,~,uz_s] = cyclic_map_seam(lat_g,lon_g,uzmo_all);
[~,~,oz_s] = cyclic_map_seam(lat_g,lon_g,ozmo_all);

%% log10-trans
cz_r = log10(cz_s+1e-4);
mz_r = log10(mz_s+1e-4);
nz_r = log10(nz_s+1e-4);
gz_r = log10(gz_s+1e-4);
iz_r = log10(iz_s+1e-4);
uz_r = log10(uz_s+1e-4);
oz_r = log10(oz_s+1e-4);

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

%% log10 Map
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(oz_r))
crameri('batlow')
%cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
cb = colorbar;%('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso (log_1_0 mgC m^-^2)')

%B - CAN
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(cz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')

%E -
% subplot('Position',[0.5 0.75 0.44 0.25])

%F - CMCC
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(mz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')

%G - GFDL
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(gz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_r))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 4])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
stamp(stex)
print('-dpng',[ppath 'Map_all_hist_clim_glmm100_log10_batlow.png'])

%% -1 to 1 Map
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - obsGLMM shifted
subplot('Position',[0.285 0.75 0.475 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(oz_n))
crameri('batlow')
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
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')

%C - CNRM
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(nz_n))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')

%D - IPSL
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(iz_n))
crameri('batlow')
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
crameri('batlow')
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
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')

%H - UK
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_s,lon_s,(uz_n))
crameri('batlow')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1])
set(gcf,'renderer','painters')
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
stamp(stex)
print('-dpng',[ppath 'Map_all_hist_clim_glmm100_neg11_scale_log10.png'])

%% Bias log10 trans Map
figure(3)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_cr)
cmocean('balance')
caxis([-2 2])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_mr)
cmocean('balance')
caxis([-2 2])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nr)
cmocean('balance')
caxis([-2 2])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gr)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso bias (log_1_0 mgC m^-^2)')
caxis([-2 2])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ir)
cmocean('balance')
caxis([-2 2])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ur)
cmocean('balance')
caxis([-2 2])
text(0.2,1.65,'UK','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%stamp(stex)
print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm100_log10.png'])

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
xlabel(cb,'bias log_1_0(zmeso) scaled [-1 1]')
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
print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm100_neg11_scale_log10.png'])

%% log10 bias with binned colorbar
cmap = cmocean('balance',7);

figure(5)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_cr)
colormap(cmap)
caxis([-1.5 1.5])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_mr)
colormap(cmap)
caxis([-1.5 1.5])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_nr)
colormap(cmap)
caxis([-1.5 1.5])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_gr)
colormap(cmap)
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'bias zmeso (log_1_0 mgC m^-^2)')
caxis([-1.5 1.5])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ir)
colormap(cmap)
caxis([-1.5 1.5])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,diff_ur)
colormap(cmap)
caxis([-1.5 1.5])
text(0.2,1.65,'UK','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%stamp(stex)
print('-dpng',[ppath 'Map_bias_all_hist_clim_glmm100_log10_bins.png'])
