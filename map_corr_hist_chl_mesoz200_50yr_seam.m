% CMIP6 output
% mesoz 200m integrations
% Map ESM future diffs 50-yr means
% And correlation mesoz-chl Hist 50-yr

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% chl corr hist
load('cmip6_hist_ssp585_corr_50yr_zmeso200_schl_sst_same_orientation',...
    'cc_s','mc_s','nc_s','gc_s','ic_s','uc_s','lat_corr','lon_corr')

% delta same orientation
load('cmip6_hist_ssp585_space_means_50yr_zmeso200_schl_sst_same_orientation.mat');

diff_ct = FCsst50 - HCsst50;
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

%% Corr
%f5 = figure('Units','inches','Position',[1 3 7.5 10]);

figure(5)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cc_s)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc_s)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CMCC','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc_s)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc_s)
cmocean('balance')
cb = colorbar('Position',[0.8 0.3 0.03 0.5],'orientation','vertical');
xlabel(cb,'zmeso-chl corr')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic_s)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc_s)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UK','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_chl_zmeso_corr_Hist.png'])

%% Corr
f5 = figure('Units','inches','Position',[1 3 6 4]);
%1 - Hist cmcc
subplot('Position',[0.01 0.45 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist cnrm
subplot('Position',[0.01 0.1 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);


% Middle
%3
subplot('Position',[0.32 0.45 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc_s)
cmocean('balance')
caxis([-1 1])
text(0,2.2,'zmeso-chl corr','HorizontalAlignment','center','FontWeight','bold')
text(-1.95,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4
subplot('Position',[0.32 0.1 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc_s)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% Right
% 5
subplot('Position',[0.63 0.45 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic_s)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6
subplot('Position',[0.63 0.1 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc_s)
cmocean('balance')
caxis([-1 1])
text(-1.95,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.32 0.075 0.3 0.02],'orientation','horizontal');

print('-dpng',[ppath 'Map_chl_zmeso_corr_Hist_v2.png'])

