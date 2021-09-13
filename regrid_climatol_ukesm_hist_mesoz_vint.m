% Regrid CMIP6 data to be same as
% COPEPOD seasonal 
% Put units into total zooplankton carbon data (mg C m-2)

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Hist
load([hpath 'ukesm_isimip_hist_zmeso_vint_climatol_1950_2014.mat']);
load([hpath 'ukesm_isimip_depth.mat']);

[LAT,LON] = meshgrid(lat,lon);

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g

uzD = double(zmc_DJF * 12.01 * 1e3);
uzJ = double(zmc_JJA * 12.01 * 1e3);
uzM = double(zmc_MAM * 12.01 * 1e3);
uzS = double(zmc_SON * 12.01 * 1e3);
uzA = double(zmc_all * 12.01 * 1e3);

%% 1 degree
lat1 = -80:90;
lon1 = -180:180;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% 
zD = griddata(LAT,LON,uzD,lat_g,lon_g);
zM = griddata(LAT,LON,uzM,lat_g,lon_g);
zJ = griddata(LAT,LON,uzJ,lat_g,lon_g);
zS = griddata(LAT,LON,uzS,lat_g,lon_g);
zall = griddata(LAT,LON,uzA,lat_g,lon_g);

%%
figure
axesm ('Robinson','MapLatLimit',[-90 90],'MapLonLimit',[30 390],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(zall))
colormap('jet')
colorbar
caxis([2 4])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('All 1^o')
%print('-dpng','gfdl_zoo_carbon_mass_all.png')

%% save
units = 'Total Carbon Mass mgC m-2';
save([hpath 'ukesm_hist_zmeso_vint_climatol_1950_2014_gridded.mat'],...
    'lat_g','lon_g',...
    'units','zD','zM','zJ','zS','zall');






