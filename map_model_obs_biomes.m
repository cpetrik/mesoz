% Maps of ESM and obs biomes
% Hist and SSP585

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/new_Ryan_chl/';
cpath1='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';
cpath2='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/ssp585/';

load([fpath 'CanESM_historical_biomes.mat'],'hcbiomes');
load([cpath1 'CMCC_historical_biomes.mat'],'hmbiomes');
load([fpath 'CNRM_historical_biomes.mat'],'hnbiomes');
load([fpath 'GFDL_historical_biomes.mat'],'hgbiomes');
load([fpath 'IPSL_historical_biomes.mat'],'hibiomes');
load([fpath 'UKESM_historical_biomes.mat'],'hubiomes');

load([fpath 'CanESM_ssp585_biomes.mat'],'scbiomes');
load([cpath2 'CMCC_ssp585_biomes.mat'],'smbiomes');
load([fpath 'CNRM_ssp585_biomes.mat'],'snbiomes');
load([fpath 'GFDL_ssp585_biomes.mat'],'sgbiomes');
load([fpath 'IPSL_ssp585_biomes.mat'],'sibiomes');
load([fpath 'UKESM_ssp585_biomes.mat'],'lat','lon','subiomes');

[lat_g,lon_g] = meshgrid(lat,lon);
clear lat lon

%%
load([bpath 'data_biomes_MODISAqua_x1.mat']);
[lat_c,lon_c] = meshgrid(lat,lon);

%% Maps

clatlim=[-90 90];
clonlim=[-280 80];
load coastlines;

cmv = colormap(viridis(3));
close all

%% Get rid of seam
lon_s = lon_g;
lat_s = lat_g;
hcbiomes_s = hcbiomes;
hmbiomes_s = hmbiomes;
hnbiomes_s = hnbiomes;
hgbiomes_s = hgbiomes;
hibiomes_s = hibiomes;
hubiomes_s = hubiomes;
scbiomes_s = scbiomes;
smbiomes_s = smbiomes;
snbiomes_s = snbiomes;
sgbiomes_s = sgbiomes;
sibiomes_s = sibiomes;
subiomes_s = subiomes;

lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
hcbiomes_s(361,:) = hcbiomes(360,:);
hmbiomes_s(361,:) = hmbiomes(360,:);
hnbiomes_s(361,:) = hnbiomes(360,:);
hgbiomes_s(361,:) = hgbiomes(360,:);
hibiomes_s(361,:) = hibiomes(360,:);
hubiomes_s(361,:) = hubiomes(360,:);
scbiomes_s(361,:) = scbiomes(360,:);
smbiomes_s(361,:) = smbiomes(360,:);
snbiomes_s(361,:) = snbiomes(360,:);
sgbiomes_s(361,:) = sgbiomes(360,:);
sibiomes_s(361,:) = sibiomes(360,:);
subiomes_s(361,:) = subiomes(360,:);

lon_o = lon_c;
lat_o = lat_c;
biomes_o = biomes;
lon_o(361,:) = lon_c(360,:)+1;
lat_o(361,:) = lat_c(360,:);
biomes_o(361,:) = biomes(360,:);

%% Hist
close all
figure(1)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hcbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hnbiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hgbiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hibiomes_s)
colormap(cmv)
caxis([1 3])
colorbar('Position',[0.8 0.3 0.025 0.4],'Ticks',[1:3],...
    'TickLabels',{'LC','HCSS','HCPS'});
text(0,1.75,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hubiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_o,lon_o,biomes_o)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'obsGLMM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_obs_ESM_biomes.png'])

%% SSP5858
figure(2)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,scbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'CAN','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,snbiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'CNRM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,sgbiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'GFDL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,sibiomes_s)
colormap(cmv)
caxis([1 3])
colorbar('Position',[0.8 0.3 0.025 0.4],'Ticks',[1:3],...
    'TickLabels',{'LC','HCSS','HCPS'});
text(0,1.75,'IPSL','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,subiomes_s)
colormap(cmv)
%caxis([0 2.25e-3])
text(0,1.75,'UKESM','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_all_ssp585_ESM_biomes.png'])

%% figure info
f1 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist cmcc
subplot('Position',[0.025 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hmbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'Historical','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist cnrm
subplot('Position',[0.025 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hnbiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist gfdl
subplot('Position',[0.025 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hgbiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Hist IPSL
subplot('Position',[0.025 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hibiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Hist uk
subplot('Position',[0.025 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hubiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6 - Hist obs
subplot('Position',[0.025 0.0 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_o,lon_o,biomes_o)
colormap(cmv)
%caxis([1 3])
text(-1.8,1.75,'obsGLMM')%,'HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% E - Fore CMCC
subplot('Position',[0.43 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,smbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'SSP5-8.5','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%F - Fore CNRM
subplot('Position',[0.43 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,snbiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%G - Fore GFDL
subplot('Position',[0.43 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,sgbiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 3])
colorbar('Position',[0.8 0.45 0.025 0.25],'Ticks',[1:3],...
    'TickLabels',{'LC','HCSS','HCPS'});

%H - Fore IPSL
subplot('Position',[0.43 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,sibiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Fore uk
subplot('Position',[0.43 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,subiomes_s)
colormap(cmv)
%caxis([1 3])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_all_hist_ssp585_ESM_biomes_cmcc.png'])

%% CAN only
figure
%1 - Hist can
subplot('Position',[0.025 0.5 0.4 0.35])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,hcbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'Historical','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% E - Fore CAN
subplot('Position',[0.43 0.5 0.4 0.35])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,scbiomes_s)
colormap(cmv)
%caxis([1 3])
text(0,1.75,'SSP5-8.5','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% colorbar('Ticks',[1:3],...
%     'TickLabels',{'LC','HCSS','HCPS'});
colorbar('Position',[0.85 0.5 0.025 0.35],'Ticks',[1:3],...
    'TickLabels',{'LC','HCSS','HCPS'});
print('-dpng',[ppath 'Map_CAN_hist_ssp585_ESM_biomes.png'])
