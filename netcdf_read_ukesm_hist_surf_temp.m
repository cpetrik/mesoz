% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server
% Surface temp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%%
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_tos_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
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
runs = find(yr>1950 & yr<=2015);

% 
sst = tos(:,:,runs);
sst=fliplr(sst);

%%
save([fpath 'ukesm_isimip_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr','runs','time','long_name','standard_name','lat','lon',...
    'units');





