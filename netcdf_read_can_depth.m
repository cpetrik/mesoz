% Read CMIP6 netcdfs
% CAN

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';

%% 
ncdisp([fpath 'deptho_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'deptho_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 291
%lon: 360
%lev: 75
%NaNs = 1.000000020040877e+20

%%
deptho(deptho(:)>1e20) = NaN;
deptho = double(deptho);

%%
save([fpath 'can_depth.mat'],'deptho','latitude','longitude')

