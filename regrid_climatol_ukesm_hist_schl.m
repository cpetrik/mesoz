% Regrid CMIP6 data to be same as
% COPEPOD seasonal 
% Put units into total zooplankton carbon data (mg C m-2)

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% Hist
load([hpath 'ukesm_isimip_hist_chl_sst_climatol_1950_2014.mat']);
[LAT,LON] = meshgrid(lat,lon);

%% Units
csD = double(sst_DJF);
csJ = double(sst_JJA);
csM = double(sst_MAM);
csS = double(sst_SON);
csA = double(sst_all);

ccD = double(chl_DJF);
ccJ = double(chl_JJA);
ccM = double(chl_MAM);
ccS = double(chl_SON);
ccA = double(chl_all);

ccD(ccD(:)<0) = 0;
ccJ(ccJ(:)<0) = 0;
ccM(ccM(:)<0) = 0;
ccS(ccS(:)<0) = 0;
ccA(ccA(:)<0) = 0;

%% chl is flipped and LAT, LON
ccD = fliplr(ccD);
ccJ = fliplr(ccJ);
ccM = fliplr(ccM);
ccS = fliplr(ccS);
ccA = fliplr(ccA);

LAT = fliplr(LAT);
LON = fliplr(LON);

%%
figure
pcolor(LAT)
colorbar

figure
pcolor(LON)
colorbar

figure
pcolor(csA)
colorbar

figure
pcolor(ccA)
colorbar

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% Use interp2 for data that are already on a regular grid
sD = interp2(LAT,LON,csD,lat_g,lon_g);
sM = interp2(LAT,LON,csM,lat_g,lon_g);
sJ = interp2(LAT,LON,csJ,lat_g,lon_g);
sS = interp2(LAT,LON,csS,lat_g,lon_g);
sall = interp2(LAT,LON,csA,lat_g,lon_g);

cD = interp2(LAT,LON,ccD,lat_g,lon_g);
cM = interp2(LAT,LON,ccM,lat_g,lon_g);
cJ = interp2(LAT,LON,ccJ,lat_g,lon_g);
cS = interp2(LAT,LON,ccS,lat_g,lon_g);
call = interp2(LAT,LON,ccA,lat_g,lon_g);

%%
figure
axesm ('Robinson','MapLatLimit',[-90 90],'MapLonLimit',[30 390],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,sall)
colormap('jet')
colorbar
caxis([0 30])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('All 1^o')
print('-dpng','ukesm_sst_all.png')

%%
figure
axesm ('Robinson','MapLatLimit',[-90 90],'MapLonLimit',[30 390],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(call))
colormap('jet')
colorbar
caxis([-7.5 -5.5])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('All 1^o')
print('-dpng','ukesm_chl_all.png')

%% save
save([hpath 'ukesm_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'lat_g','lon_g',...
    'units','cD','cM','cJ','cS','call','sD','sM','sJ','sS','sall');






