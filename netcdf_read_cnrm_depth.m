% Read CMIP6 netcdfs
% CNRM

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% 
ncdisp([fpath 'deptho_Ofx_CNRM-ESM2-1_historical_r2i1p1f2_gn.nc'])

%%
ncid = netcdf.open([fpath 'deptho_Ofx_CNRM-ESM2-1_historical_r2i1p1f2_gn.nc'],'NC_NOWRITE');

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

%%
depth = deptho;
depth(depth(:)>1e20) = NaN;
wet = find(~isnan(depth(:,:,1)));
depv = reshape(depth,362*294,75);
depv = depv(wet,:);
dep = NaN*ones(size(wet));

%%
for z=1:length(dep)
    bot = find(isnan(depv(z,:)));
    dep(z) = depv(z,bot(1)-1);
end

deptho = NaN*ones(362,294);
deptho(wet) = dep;

%%
save([fpath 'cnrm_depth.mat'],'deptho')

