% CMIP6 output 
% mesoz 200m integrations, surf chl, stt

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat']);

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
load([mpath 'cmcc_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat'])
 
%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat']);

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat'])

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat']);

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat']);

%%
clatlim=[-90 90];
clonlim=[-280 80];

vpath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([vpath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon');
[lat_g,lon_g] = meshgrid(lat,lon);

load coastlines;

%% chl
figure
pcolor(cccorr)
shading flat
title('CAN')

figure
pcolor(mccorr)
shading flat
title('CMCC')

figure
pcolor(nccorr)
shading flat
title('CNRM')

figure
pcolor(iccorr)
shading flat
title('IP')

figure
pcolor(gccorr)
shading flat
title('GF')

figure
pcolor(uccorr)
shading flat
title('UK')

%% sst
close all
figure
pcolor(ctcorr)
shading flat
title('CAN')

figure
pcolor(mtcorr)
shading flat
title('CMCC')

figure
pcolor(ntcorr)
shading flat
title('CNRM')

figure
pcolor(itcorr)
shading flat
title('IP')

figure
pcolor(gtcorr)
shading flat
title('GF')

figure
pcolor(utcorr)
shading flat
title('UK')

figure
pcolor(lat_g)
shading flat
title('lat')

%% flip as needed
close all
gccorr = fliplr(gccorr);
gtcorr = fliplr(gtcorr);
gccorr2 = fliplr(gccorr2);
gtcorr2 = fliplr(gtcorr2);

%% Fix seam
lon_s = lon_g;
lat_s = lat_g;

cc_s = cccorr;
mc_s = mccorr;
nc_s = nccorr;
gc_s = gccorr;
ic_s = iccorr;
uc_s = uccorr;

cc2_s = cccorr2;
mc2_s = mccorr2;
nc2_s = nccorr2;
gc2_s = gccorr2;
ic2_s = iccorr2;
uc2_s = uccorr2;

lon_s(361,:) = lon_s(360,:)+1;
lat_s(361,:) = lat_s(360,:);
cc_s(361,:) = cccorr(360,:)+eps;
mc_s(361,:) = mccorr(360,:)+eps;
nc_s(361,:) = nccorr(360,:)+eps;
gc_s(361,:) = gccorr(360,:)+eps;
ic_s(361,:) = iccorr(360,:)+eps;
uc_s(361,:) = uccorr(360,:)+eps;

cc2_s(361,:) = cccorr2(360,:)+eps;
mc2_s(361,:) = mccorr2(360,:)+eps;
nc2_s(361,:) = nccorr2(360,:)+eps;
gc2_s(361,:) = gccorr2(360,:)+eps;
ic2_s(361,:) = iccorr2(360,:)+eps;
uc2_s(361,:) = uccorr2(360,:)+eps;

ct_s = ctcorr;
mt_s = mtcorr;
nt_s = ntcorr;
gt_s = gtcorr;
it_s = itcorr;
ut_s = utcorr;

ct2_s = ctcorr2;
mt2_s = mtcorr2;
nt2_s = ntcorr2;
gt2_s = gtcorr2;
it2_s = itcorr2;
ut2_s = utcorr2;

ct_s(361,:) = ctcorr(360,:)+eps;
mt_s(361,:) = mtcorr(360,:)+eps;
nt_s(361,:) = ntcorr(360,:)+eps;
gt_s(361,:) = gtcorr(360,:)+eps;
it_s(361,:) = itcorr(360,:)+eps;
ut_s(361,:) = utcorr(360,:)+eps;

ct2_s(361,:) = ctcorr2(360,:)+eps;
mt2_s(361,:) = mtcorr2(360,:)+eps;
nt2_s(361,:) = ntcorr2(360,:)+eps;
gt2_s(361,:) = gtcorr2(360,:)+eps;
it2_s(361,:) = itcorr2(360,:)+eps;
ut2_s(361,:) = utcorr2(360,:)+eps;

%% save same orient
lat_corr = lat_s;
lon_corr = lon_s;
save('cmip6_hist_ssp585_corr_50yr_zmeso200_schl_sst_same_orientation',...
    'cc_s','mc_s','nc_s','gc_s','ic_s','uc_s',...
    'cc2_s','mc2_s','nc2_s','gc2_s','ic2_s','uc2_s',...
    'ct_s','mt_s','nt_s','gt_s','it_s','ut_s',...
    'ct2_s','mt2_s','nt2_s','gt2_s','it2_s','ut2_s','lat_corr','lon_corr')

%% test function
[lat2,lon2,gvar2] = cyclic_map_seam(lat_g,lon_g,gtcorr);

figure
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gt_s)
cmocean('balance')
caxis([-1 1])

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat2,lon2,gvar2)
cmocean('balance')
caxis([-1 1])

%% Just Hist
close all
figure(10)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cc_s)
cmocean('balance')
caxis([-1 1])
%text(0,1.75,'Historic chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%colorbar

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_all_hist_chl_mesoz_corr.png'])

%% Map Chl
close all
% figure info
f1 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist cmcc
subplot('Position',[0.025 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'Historic chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist cnrm
subplot('Position',[0.025 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist gfdl
subplot('Position',[0.025 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Hist IPSL
subplot('Position',[0.025 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Hist uk
subplot('Position',[0.025 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6 - Hist obs
% subplot('Position',[0.025 0.0 0.4 0.165])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(lat_o,lon_o,biomes_o)
% cmocean('balance')
% caxis([-1 1])
% text(-1.75,1.75,'obs','HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 7 - Fore CMCC
subplot('Position',[0.43 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mc2_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'SSP5-8.5 chl','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Fore CNRM
subplot('Position',[0.43 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nc2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - Fore GFDL
subplot('Position',[0.43 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gc2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.8 0.45 0.025 0.25]);

%10 - Fore IPSL
subplot('Position',[0.43 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ic2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%11 - Fore uk
subplot('Position',[0.43 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uc2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_all_hist_ssp585_chl_mesoz_corr.png'])

%% figure info SST
f2 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist cmcc
subplot('Position',[0.025 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mt_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'Historic SST','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist cnrm
subplot('Position',[0.025 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nt_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist gfdl
subplot('Position',[0.025 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gt_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Hist IPSL
subplot('Position',[0.025 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,it_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - Hist uk
subplot('Position',[0.025 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ut_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6 - Hist obs
% subplot('Position',[0.025 0.0 0.4 0.165])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(lat_o,lon_o,biomes_o)
% cmocean('balance')
% caxis([-1 1])
% text(-1.75,1.75,'obs','HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 7 - Fore CMCC
subplot('Position',[0.43 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mt2_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'SSP5-8.5 SST','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CMCC','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Fore CNRM
subplot('Position',[0.43 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nt2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - Fore GFDL
subplot('Position',[0.43 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gt2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.8 0.45 0.025 0.25]);

%10 - Fore IPSL
subplot('Position',[0.43 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,it2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%11 - Fore uk
subplot('Position',[0.43 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ut2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_all_hist_ssp585_sst_mesoz_corr.png'])

%% No CAN
f3 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist cnrm
subplot('Position',[0.025 0.73 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nt_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'Historic SST','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%1 - Hist gfdl
subplot('Position',[0.025 0.49 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gt_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist IPSL
subplot('Position',[0.025 0.25 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,it_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - Hist uk
subplot('Position',[0.025 0.01 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ut_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);


%5 - Fore CNRM
subplot('Position',[0.43 0.73 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nt2_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'SSP5-8.5 SST','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6 - Fore GFDL
subplot('Position',[0.43 0.49 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gt2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.85 0.35 0.025 0.3]);

%7 - Fore IPSL
subplot('Position',[0.43 0.25 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,it2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - Fore uk
subplot('Position',[0.43 0.01 0.4 0.24])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ut2_s)
cmocean('balance')
caxis([-1 1])
text(-1.75,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_all_hist_ssp585_sst_mesoz_corr_noCAN.png'])

%% Just CAN for supp
% figure info
figure(5) %('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist can chl
subplot('Position',[0.025 0.5 0.4 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cc_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'Historic chl','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist can sst
subplot('Position',[0.025 0.1 0.4 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ct_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'Historic SST','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 7 - Fore CAN chl
subplot('Position',[0.43 0.5 0.4 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cc2_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'SSP5-8.5 chl','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.85 0.55 0.025 0.3]);

%8 - Fore can sst
subplot('Position',[0.43 0.1 0.4 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ct2_s)
cmocean('balance')
caxis([-1 1])
text(0,1.75,'SSP5-8.5 SST','HorizontalAlignment','center','FontWeight','bold')
%text(-1.75,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar('Position',[0.85 0.15 0.025 0.3]);

print('-dpng',[ppath 'Map_CAN_hist_ssp585_chl_sst_mesoz_corr.png'])
