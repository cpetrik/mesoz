% Regrid CMIP6 data to be same as
% COPEPOD seasonal 
% Put units into total zooplankton carbon data (mg C m-2)

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';


%% Hist
load([hpath 'cnrm_hist_zmeso_vint_climatol_1950_2014.mat']);
load([hpath 'cnrm_depth.mat'],'deptho');

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g

czD = double(zmc_DJF) * 12.01 * 1e3;
czJ = double(zmc_JJA) * 12.01 * 1e3;
czM = double(zmc_MAM) * 12.01 * 1e3;
czS = double(zmc_SON) * 12.01 * 1e3;
czA = double(zmc_all) * 12.01 * 1e3;

%% 1 degree
lat1 = -80:90;
lon1 = -180:180;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% 
LAT = lat;
LON = lon;

zD = griddata(LAT,LON,czD,lat_g,lon_g);
zM = griddata(LAT,LON,czM,lat_g,lon_g);
zJ = griddata(LAT,LON,czJ,lat_g,lon_g);
zS = griddata(LAT,LON,czS,lat_g,lon_g);
zall = griddata(LAT,LON,czA,lat_g,lon_g);

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
save([hpath 'cnrm_hist_zmeso_vint_climatol_1950_2014_gridded.mat'],...
    'lat_g','lon_g',...
    'units','zD','zM','zJ','zS','zall');






