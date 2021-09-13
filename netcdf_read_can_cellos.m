% Read CMIP6 netcdfs
% CAN

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';

%% Cell thickness
ncdisp([fpath 'thkcello_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'thkcello_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 291
%lon: 360
%lev: 45
%NaNs = 1.000000020040877e+20

%%
thkcello(thkcello(:)>1e20) = NaN;

%% Cell area
ncdisp([fpath 'areacello_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'areacello_Ofx_CanESM5-CanOE_historical_r1i1p2f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 291
%lon: 360
%NaNs = 1.000000020040877e+20

areacello(areacello(:)>1e20) = NaN;

%% Cell volume
ncdisp([fpath 'volcello_Ofx_CanESM5-CanOE_esm-hist_r1i1p2f1_gn.nc'])
ncdisp([fpath 'volcello_Ofx_CanESM5-CanOE_omip1_r1i1p2f1_gn.nc'])
ncdisp([fpath 'volcello_Ofx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'volcello_Ofx_CanESM5-CanOE_omip1_r1i1p2f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 291
%lon: 360
%lev: 45
%NaNs = 1.000000020040877e+20

volcello(volcello(:)>1e20) = NaN;

%% test
areac = repmat(areacello,1,1,45);
tv = thkcello(~isnan(thkcello(:)));
av = areac(~isnan(thkcello(:)));
vv = volcello(~isnan(volcello(:)));

thk = vv ./ av;
eq = (thk==tv);

%%
load([fpath 'can_depth.mat'],'deptho')

save([fpath 'can_grid.mat'],'deptho','latitude','longitude',...
    'thkcello','areacello','lev','lev_bnds')

