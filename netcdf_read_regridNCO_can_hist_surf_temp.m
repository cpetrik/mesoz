% Read CMIP6 CanESM5 Historic netcdfs
% surface temp
% Regridded with NCO

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

%% sst
ncdisp([fpath 'CanESM5-CanOE_r1i1p2f1_historical_tos_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
% comment     = 'Temperature of upper boundary of the liquid ocean, including temperatures below sea-ice and floating ice shelves.'
units         = 'degC';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x291x1980
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'CanESM5-CanOE_r1i1p2f1_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test=tos(:,:,1900);
%%
test(test>1e19)=NaN;

pcolor(test)

min(test(:))
max(test(:))

