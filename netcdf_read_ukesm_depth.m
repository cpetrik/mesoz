% Read CMIP6 netcdfs
% UK-ESM Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Meso Zoop zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_picontrol_deptho_onedeg_global_fx.nc'])

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_picontrol_deptho_onedeg_global_fx.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 360
%lon: 180
%NaNs = 1.000000020040877e+20

%%
save([fpath 'ukesm_depth.mat'],'deptho')