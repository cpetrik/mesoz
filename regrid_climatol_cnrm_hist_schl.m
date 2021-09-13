% Regrid CMIP6 data to be same as
% COPEPOD seasonal 
% Put units into total zooplankton carbon data (mg C m-2)

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% Hist
load([hpath 'cnrm_hist_chl_sst_climatol_1950_2014.mat']);

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

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% input data not regular (tripolar?)
LAT = lat;
LON = lon;

sD = griddata(LAT,LON,csD,lat_g,lon_g);
sM = griddata(LAT,LON,csM,lat_g,lon_g);
sJ = griddata(LAT,LON,csJ,lat_g,lon_g);
sS = griddata(LAT,LON,csS,lat_g,lon_g);
sall = griddata(LAT,LON,csA,lat_g,lon_g);

cD = griddata(LAT,LON,ccD,lat_g,lon_g);
cM = griddata(LAT,LON,ccM,lat_g,lon_g);
cJ = griddata(LAT,LON,ccJ,lat_g,lon_g);
cS = griddata(LAT,LON,ccS,lat_g,lon_g);
call = griddata(LAT,LON,ccA,lat_g,lon_g);

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
print('-dpng','cnrm_sst_all.png')

figure
axesm ('Robinson','MapLatLimit',[-90 90],'MapLonLimit',[30 390],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(call))
colormap('jet')
colorbar
caxis([-4.5 -3])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('All 1^o')
print('-dpng','cnrm_chl_all.png')

%% save
save([hpath 'cnrm_hist_chl_sst_climatol_1950_2014_gridded.mat'],...
    'lat_g','lon_g',...
    'units','cD','cM','cJ','cS','call','sD','sM','sJ','sS','sall');







