% Regrid CMIP6 data to regular 1 degree grid

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';

%% Hist
load([hpath 'can_ssp585_chl_sst_climatol_2051_2100.mat']);

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

%%
figure
pcolor(latitude)
colorbar

figure
pcolor(longitude)
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

%% longitude is shifted 0-360 and grid is not regulat (tripolar?)
test = longitude; %-360;
id=find(test > 180);
test(id)=test(id)-360;
LON = double(test);

LAT = latitude;
%LON = longitude;

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
%print('-dpng','can_sst_all.png')

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
%print('-dpng','can_chl_all.png')

%% save
save([hpath 'can_ssp585_chl_sst_climatol_2051_2100_gridded.mat'],...
    'lat_g','lon_g',...
    'units','cD','cM','cJ','cS','call','sD','sM','sJ','sS','sall');






