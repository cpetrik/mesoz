% Read CMIP6 netcdfs
% CNRM

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% Cell area
ncdisp([fpath 'areacello_Ofx_CNRM-ESM2-1_historical_r1i1p1f2_gn.nc'])

%%
ncid = netcdf.open([fpath 'areacello_Ofx_CNRM-ESM2-1_historical_r1i1p1f2_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 294
%lon: 362
%lev: 75
%NaNs = 1.000000020040877e+20

areacello(areacello(:)>1e20) = NaN;

%% Cell volume?
ncdisp([fpath 'volo_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc'])

%just a time series of total volume

%% Cell thickness
%Varies over time, just get directly in integration

%%
load([fpath 'cnrm_depth.mat'],'deptho')

save([fpath 'cnrm_grid.mat'],'deptho','lat','lon',...
    'areacello')

