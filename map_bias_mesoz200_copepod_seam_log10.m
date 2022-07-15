% CMIP6 output
% mesoz 200m integrations
% Map ESM bias raw and scaled
% log-10 trans

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

Cdir = '/Volumes/MIP/Fish-MIP/CMIP6/';
load([Cdir 'GFDL/gridspec_gfdl_cmip6.mat'],'deptho','LAT','LON','lmask');

%% load copepod gridded
load('copepod-2012_cmass_all_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

lat_c = lat_g; 
lon_c = lon_g;
clear lat_g lon_g

%%
load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat');

stex = 'log10';

%% COPEPOD units from mg m-3 to mg m-2 (integrate top 200m)
dep200 = min(200,deptho);
dep200(isnan(deptho(:))) = nan;
zoo_200 = zoo_g .* dep200;

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
czmo_all = czmo_all * 12.01 * 1e3;
mzmo_all = mzmo_all * 12.01 * 1e3;
nzmo_all = nzmo_all * 12.01 * 1e3;
gzmo_all = gzmo_all * 12.01 * 1e3;
izmo_all = izmo_all * 12.01 * 1e3;
uzmo_all = uzmo_all * 12.01 * 1e3;
%obsCOPE mg C/m2

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
[~,~,oz_s] = cyclic_map_seam(lat_g,lon_g,zoo_200);

%% log10-trans
cz_r = log10(cz_s+1e-4);
mz_r = log10(mz_s+1e-4);
nz_r = log10(nz_s+1e-4);
gz_r = log10(gz_s+1e-4);
iz_r = log10(iz_s+1e-4);
uz_r = log10(uz_s+1e-4);
oz_r = log10(oz_s+1e-4);

%% Diffs (Bias)
diff_cr = cz_r - oz_r;
diff_mr = mz_r - oz_r;
diff_nr = nz_r - oz_r;
diff_gr = gz_r - oz_r;
diff_ir = iz_r - oz_r;
diff_ur = uz_r - oz_r;

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
print('-dpng',[ppath 'Map_bias_all_hist_clim_copepod_log10.png'])

