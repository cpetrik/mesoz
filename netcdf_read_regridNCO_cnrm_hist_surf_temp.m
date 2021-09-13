% Read CMIP6 netcdfs
% CNRM-ESM2-1 Hist
% Surface temp
% Regridded with NCO

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';

%%
ncdisp([fpath 'CNRM-ESM2-1_r1i1p1f2_historical_tos_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
% cell_measures = 'area: areacello'
%Size:       362x294x1980
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'CNRM-ESM2-1_r1i1p1f2_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test=tos(:,:,1900);
test(test>1e19)=NaN;

pcolor(test)

min(test(:)) %TEMP SPANS -14 TO 17
