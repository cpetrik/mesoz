% Read CMIP6 netcdfs
% UK-ESM SSP 126
% NetCDF on the DKRZ ISIMIP server
% Surface temp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp126/';

%%
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_ssp126_tos_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_ssp126_tos_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
tos(tos>=1e20) = nan;

%% Time
yr = ((time+1)/12)+1601;

% 
sst = tos;
sst=fliplr(sst);

%%
save([fpath 'ukesm_isimip_ssp126_sst_monthly_2015_2100.mat'],'sst',...
    'yr','time','long_name','standard_name','lat','lon',...
    'units');





