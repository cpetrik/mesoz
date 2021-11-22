% Read Seawifs chl biomes for use with Stromberg 

clear all
close all

cpath='/Volumes/MIP/Fish-MIP/CMIP6/biome_masks/';
fpath='/Volumes/MIP/Obs_Data/Chl/';

%% biomes
ncdisp([cpath 'SeaWiFS_based_biomes_x1.nc'])

%%
%biomes
% Size:       360x180
% Dimensions: lon,lat
% Datatype:   double
% _FillValue            = NaN
key                   = '1-LC, 2-HCSS, 3-HCPS';
chlorophyll_threshold = '0.125 mg Chl m-3';

%%
ncid = netcdf.open([cpath 'SeaWiFS_based_biomes_x1.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
sbiomes = biomes;

%%
save([cpath 'SeaWiFS_based_biomes_x1.mat'],'sbiomes',...
    'lat','lon','key','chlorophyll_threshold');
save([fpath 'SeaWiFS_based_biomes_x1.mat'],'sbiomes',...
    'lat','lon','key','chlorophyll_threshold');





