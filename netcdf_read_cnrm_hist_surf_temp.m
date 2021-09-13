% Read CMIP6 netcdfs
% CNRM-ESM2-1 Hist
% Surface temp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';

%%
ncdisp([fpath 'tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
% cell_measures = 'area: areacello'
%Size:       362x294x1980
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'tos_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time+1)/365)+1850;
runs = find(yr>1950 & yr<=2015);
sst = tos(:,:,runs);
sst(sst >= 1.00e+18) = NaN;

%%
save([fpath 'cnrm_hist_sst_monthly_1950_2014.mat'],'sst',...
    'long_name','standard_name','units',...
    'lat','lon','time','time_bounds','yr','runs');






