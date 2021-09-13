% Read Biomes def by Stock et al. 2014
% Calc by J Luo from MODIS obs (used in GLMM) provided by Ryan

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';
ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% biomes
ncdisp([bpath 'data_biomes_MODISAqua_x1.nc'])

%%
% biomes
% Size:       360x180
% Dimensions: lon,lat
% Datatype:   double
% Attributes:
% _FillValue = NaN
% key        = '1-LC; 2-HCSS; 3-HCPS'

%%
ncid = netcdf.open([bpath 'data_biomes_MODISAqua_x1.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
figure
pcolor(biomes)
shading flat
colorbar

%%
save([bpath 'data_biomes_MODISAqua_x1.mat']);
load([bpath 'data_biomes_MODISAqua_x1.mat']);

%% Maps
[lat_g,lon_g] = meshgrid(lat,lon);

clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

cmv = colormap(viridis(3));

%%
figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,biomes)
colormap(cmv)
%caxis([1 3])
% colorbar('position',[0.5 0.25 0.025 0.5],'Ticks',[1:3],...
%          'TickLabels',{'LC','HCSS','HCPS'})
colorbar('Ticks',[1:3],'TickLabels',{'LC','HCSS','HCPS'})
title('MODIS')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%print('-dpng',[ppath 'Map_MODIS_biomes_hist.png'])

