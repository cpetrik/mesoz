% Read obs chl 
% 1997-2018 SeaWiFS and MODIS product + 
% the Johnson et al. Southern Ocean modification
% from J Luo

clear all
close all

fpath='/Volumes/MIP/Obs_Data/Chl/';

%% biomes
ncdisp([fpath 'globcolour_soblend_mc_x1.nc'])

%%
% chl
% Size:       360x180x12
% Dimensions: lon,lat,time
% Datatype:   double
% Attributes:
% _FillValue    = NaN
% regrid_method = 'bilinear'
% coordinates   = 'month'

%%
ncid = netcdf.open([fpath 'globcolour_soblend_mc_x1.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test = squeeze(chl(:,:,6));
pcolor(test)
shading flat
colorbar

%%
save([fpath 'globcolour_soblend_mc_x1.mat']);

