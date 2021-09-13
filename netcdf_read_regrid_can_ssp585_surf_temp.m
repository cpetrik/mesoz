% Read CMIP6 CanESM5 SSP 585 netcdfs
% surface temp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';

%% 
ncdisp([fpath 'tos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
% comment     = 'Temperature of upper boundary of the liquid ocean, including temperatures below sea-ice and floating ice shelves.'
units         = 'degC';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x180x1032
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'tos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time)/365)+1850;
sst = tos;

%%
save([fpath 'can_ssp585_sst_monthly_onedeg_2015_2100.mat'],'sst',...
    'long_name','standard_name','units',...
    'lat','lon','time','yr');
