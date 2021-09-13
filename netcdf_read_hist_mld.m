% Read MLD from each ESM CMIP6

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';

%% mixed layer depth defined by sigma t
ncdisp([fpath 'mlotst_Omon_GFDL-ESM4_historical_r1i1p1f1_gn_195001-196912.nc'])

%%
% mlotst
% Size:       720x576x240
% Dimensions: x,y,time
% Datatype:   single
% Attributes:
long_name     = 'Ocean Mixed Layer Thickness Defined by Sigma T';
units         = 'm';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
% standard_name = 'ocean_mixed_layer_thickness_defined_by_sigma_t'
% coordinates   = 'lat lon'
% original_name = 'mlotst'

%%
ncid = netcdf.open([fpath 'mlotst_Omon_GFDL-ESM4_historical_r1i1p1f1_gn_195001-196912.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
pcolor(squeeze(mlotst(:,:,6)))
shading flat
colorbar
caxis([0 200])

%%
%save([fpath 'gfdl_mld_hist_1950_2105.mat']);

