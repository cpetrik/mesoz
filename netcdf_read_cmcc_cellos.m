% Read CMIP6 netcdfs
% CMCC calc thkcello from volcello/areacello

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';

%% HISTORIC ---------------------------------------
%% Cell area
ncdisp([fpath 'hist/areacello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'hist/areacello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 292
%lon: 362
%NaNs = 1.000000020040877e+20
%ocnBgchem: BFM5.1

areacello(areacello(:)>1e20) = NaN;
lat_a = latitude;
lon_a = longitude;

%% Cell volume
ncdisp([fpath 'hist/volcello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'hist/volcello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.nc'],'NC_NOWRITE');

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
%lev: 50
%NaNs = 1.000000020040877e+20

volcello(volcello(:)>1e20) = NaN;
lat = latitude;
lon = longitude;

%% Check same orientation
figure
pcolor(volcello(:,:,1))

figure
pcolor(volcello(:,:,40))

figure
pcolor(areacello(:,:,1))

%% Cell layer thickness
arearep = repmat(areacello,1,1,50);

figure
pcolor(arearep(:,:,40))

thkcello = volcello ./ arearep;

%%
save([fpath 'hist/thkcello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.mat'],'lat','lon',...
    'thkcello','lev')

%% SSP 585 -------------------------------
clear areacello volcello thkcello lat lon arearep
%% Cell area
ncdisp([fpath 'ssp585/areacello_Ofx_CMCC-ESM2_ssp585_r1i1p1f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'ssp585/areacello_Ofx_CMCC-ESM2_ssp585_r1i1p1f1_gn.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 292
%lon: 362
%NaNs = 1.000000020040877e+20
%ocnBgchem: BFM5.1

areacello(areacello(:)>1e20) = NaN;
lat_a = latitude;
lon_a = longitude;

%% Cell volume
ncdisp([fpath 'ssp585/volcello_Ofx_CMCC-ESM2_ssp585_r1i1p1f1_gn.nc'])

%%
ncid = netcdf.open([fpath 'ssp585/volcello_Ofx_CMCC-ESM2_ssp585_r1i1p1f1_gn.nc'],'NC_NOWRITE');

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
%lev: 50
%NaNs = 1.000000020040877e+20

volcello(volcello(:)>1e20) = NaN;
lat = latitude;
lon = longitude;

%% Check same orientation
figure
pcolor(volcello(:,:,1))

figure
pcolor(volcello(:,:,40))

figure
pcolor(areacello(:,:,1))

%% Cell layer thickness
arearep = repmat(areacello,1,1,50);

figure
pcolor(arearep(:,:,40))

thkcello = volcello ./ arearep;

%%
save([fpath 'ssp585/thkcello_Ofx_CMCC-ESM2_ssp585_r1i1p1f1_gn.mat'],'lat','lon',...
    'thkcello','lev')
