% Regrid CMIP6 data to be same as
% COPEPOD seasonal 
% Put units into total zooplankton carbon data (mg C m-2)

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

%% Hist
load([hpath 'ipsl_ssp585_zmeso_200_climatol_2051_2100.mat']);

[LAT,LON] = meshgrid(lat,lon);

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to g(WW) m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g

izD = double(zmc_DJF) * 12.01 * 1e3;
izJ = double(zmc_JJA) * 12.01 * 1e3;
izM = double(zmc_MAM) * 12.01 * 1e3;
izS = double(zmc_SON) * 12.01 * 1e3;
izA = double(zmc_all) * 12.01 * 1e3;

%% flipped and LAT, LON

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
pcolor(izA)
colorbar

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% 

zD = interp2(LAT,LON,izD,lat_g,lon_g);
zM = interp2(LAT,LON,izM,lat_g,lon_g);
zJ = interp2(LAT,LON,izJ,lat_g,lon_g);
zS = interp2(LAT,LON,izS,lat_g,lon_g);
zall = interp2(LAT,LON,izA,lat_g,lon_g);

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
print('-dpng','ipsl_ssp585_zoo_carbon_mass_all.png')

%% save
units = 'Total Carbon Mass mgC m-2';
save([hpath 'ipsl_ssp585_zmeso_200_climatol_2051_2100_gridded.mat'],...
    'lat_g','lon_g',...
    'units','zD','zM','zJ','zS','zall');






