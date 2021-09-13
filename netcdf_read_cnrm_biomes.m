% Read CNRM biomes from J Luo
% based on MODIS obs (used in GLMM) provided by Ryan

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';
ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
bpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/';

%% biomes
ncdisp([bpath 'CNRM_historical_biomes_x1.nc'])

%%
% dimensions:
% lat = 180 ;
% lon = 360 ;

% variables:
% double biomes(lat, lon) ;
% biomes:_FillValue = NaN ;
% biomes:key = "1-LC; 2-HCSS; 3-HCPS" ;
% biomes:chlorophyll-threshold = "0.165 mg Chl m^-3" ;


%% Hist
ncid = netcdf.open([bpath 'CNRM_historical_biomes_x1.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
hnbiomes = biomes;
pcolor(hnbiomes)
shading flat
colorbar

%%
save([bpath 'CNRM_historical_biomes.mat'],'lat','lon','biomes','hnbiomes');

%% SSP585
clear lat lon biomes

ncid = netcdf.open([bpath 'CNRM_ssp585_biomes_x1.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
snbiomes = biomes;
figure
pcolor(snbiomes)
shading flat
colorbar

%%
save([bpath 'CNRM_ssp585_biomes.mat'],'lat','lon','biomes','snbiomes');

%% Maps
[lat_g,lon_g] = meshgrid(lat,lon);

clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

cmv = colormap(viridis(3));

%% 
close all
figure(2)
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,hnbiomes)
colormap(cmv)
%caxis([1 3])
colorbar('position',[0.5 0.25 0.025 0.5],'Ticks',[1:3],...
         'TickLabels',{'LC','HCSS','HCPS'})
title('CNRM Hist')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.05 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,snbiomes)
colormap(cmv)
%caxis([1 3])
title('CNRM SSP585')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_CNRM_biomes_hist_ssp585.png'])
