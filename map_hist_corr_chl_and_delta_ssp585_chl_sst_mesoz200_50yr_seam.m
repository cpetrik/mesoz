% CMIP6 output
% mesoz 200m integrations
% Map ESM future diffs 50-yr means
% And correlation mesoz-chl Hist 50-yr

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% delta same orientation
load('cmip6_hist_ssp585_space_means_50yr_zmeso200_schl_sst_same_orientation.mat');
% chl corr hist
load('cmip6_hist_ssp585_corr_50yr_zmeso200_schl_sst_same_orientation',...
    'cc_s','mc_s','nc_s','gc_s','ic_s','uc_s','lat_corr','lon_corr')

%% Diffs
diff_ct = FCsst50 - HCsst50;
diff_nt = FNsst50 - HNsst50;
diff_gt = FGsst50 - HGsst50;
diff_it = FIsst50 - HIsst50;
diff_ut = FUsst50 - HUsst50;

%% Pdiffs
pdiff_cc = 100* (FCchl - HCchl) ./ HCchl;
pdiff_mc = 100* (FMchl - HMchl) ./ HMchl;
pdiff_nc = 100* (FNchl - HNchl) ./ HNchl;
pdiff_gc = 100* (FGchl - HGchl) ./ HGchl;
pdiff_ic = 100* (FIchl - HIchl) ./ HIchl;
pdiff_uc = 100* (FUchl - HUchl) ./ HUchl;

pdiff_cz = 100* (FCzm - HCzm) ./ HCzm;
pdiff_mz = 100* (FMzm - HMzm) ./ HMzm;
pdiff_nz = 100* (FNzm - HNzm) ./ HNzm;
pdiff_gz = 100* (FGzm - HGzm) ./ HGzm;
pdiff_iz = 100* (FIzm - HIzm) ./ HIzm;
pdiff_uz = 100* (FUzm - HUzm) ./ HUzm;

%% Fix seam
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);
[~,~,pdiff_cc2] = cyclic_map_seam(lat_g,lon_g,pdiff_cc);
[~,~,pdiff_mc2] = cyclic_map_seam(lat_g,lon_g,pdiff_mc);
[~,~,pdiff_nc2] = cyclic_map_seam(lat_g,lon_g,pdiff_nc);
[~,~,pdiff_gc2] = cyclic_map_seam(lat_g,lon_g,pdiff_gc);
[~,~,pdiff_ic2] = cyclic_map_seam(lat_g,lon_g,pdiff_ic);
[~,~,pdiff_uc2] = cyclic_map_seam(lat_g,lon_g,pdiff_uc);

[~,~,pdiff_cz2] = cyclic_map_seam(lat_g,lon_g,pdiff_cz);
[~,~,pdiff_mz2] = cyclic_map_seam(lat_g,lon_g,pdiff_mz);
[~,~,pdiff_nz2] = cyclic_map_seam(lat_g,lon_g,pdiff_nz);
[~,~,pdiff_gz2] = cyclic_map_seam(lat_g,lon_g,pdiff_gz);
[~,~,pdiff_iz2] = cyclic_map_seam(lat_g,lon_g,pdiff_iz);
[~,~,pdiff_uz2] = cyclic_map_seam(lat_g,lon_g,pdiff_uz);

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

%% Corr
f5 = figure('Units','inches','Position',[1 3 7.5 10]);

%1 - Hist cmcc
subplot('Position',[0.01 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc_s)
cmocean('balance')
caxis([-1 1])
text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
text(0,2.2,'Historic corr','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%1 - Hist cnrm
subplot('Position',[0.01 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%1 - Hist gfdl
subplot('Position',[0.01 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist IPSL
subplot('Position',[0.01 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Hist uk
subplot('Position',[0.01 0.1 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.0195 0.075 0.275 0.02],'orientation','horizontal');

% Change
%1 - Chl cmcc
subplot('Position',[0.32 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mc2)
cmocean('balance')
caxis([-100 100])
text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
text(0,2.2,'% \Delta chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.95,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Chl cnrm
subplot('Position',[0.32 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nc2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Chl gfdl
subplot('Position',[0.32 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gc2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Chl IPSL
subplot('Position',[0.32 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_ic2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Chl uk
subplot('Position',[0.32 0.1 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uc2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - Zoo CMCC
subplot('Position',[0.63 0.8 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_mz2)
cmocean('balance')
caxis([-100 100])
text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
text(0,2.2,'% \Delta zmeso','HorizontalAlignment','center','FontWeight','bold')
text(-1.95,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%7 - Zoo CNRM
subplot('Position',[0.63 0.625 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_nz2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Zoo GFDL
subplot('Position',[0.63 0.45 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_gz2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - Zoo IPSL
subplot('Position',[0.63 0.275 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_iz2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%10 - Zoo uk
subplot('Position',[0.63 0.1 0.3 0.175])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_uz2)
cmocean('balance')
caxis([-100 100])
text(-1.95,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.475 0.075 0.3 0.02],'orientation','horizontal');

print('-dpng',[ppath 'Map_chl_zmeso_corr_pdiff_Hist_SSP585_50s_cmcc.png'])

%% Just CAN for Supp
figure(8)
%1 - Hist can chl
subplot('Position',[0.01 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_corr,lon_corr,cc_s)
cmocean('balance')
caxis([-1 1])
colorbar('Position',[0.01 0.45 0.3 0.025],'orientation','horizontal');
text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
text(0,1.75,'Historic corr','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.32 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cc2)
cmocean('balance')
caxis([-100 100])
text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
text(0,1.75,'% \Delta chl','HorizontalAlignment','center','FontWeight','bold')
%text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.63 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,pdiff_cz2)
cmocean('balance')
caxis([-100 100])
text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
text(0,1.75,'% \Delta zmeso','HorizontalAlignment','center','FontWeight','bold')
%text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.475 0.45 0.3 0.025],'orientation','horizontal');
print('-dpng',[ppath 'Map_chl_zmeso_corr_pdiff_Hist_SSP585_50s_CANonly.png'])
