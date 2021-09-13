% Read CMIP6 netcdfs
% CNRM-ESM2-1 SSP 126
% Surface temp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp126/';

%%
ncdisp([fpath 'tos_Omon_CNRM-ESM2-1_ssp126_r1i1p1f2_gn_201501-210012_onedeg.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
% cell_measures = 'area: areacello'
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'tos_Omon_CNRM-ESM2-1_ssp126_r1i1p1f2_gn_201501-210012_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time)/365)+1850;
sst = tos;

%%
save([fpath 'cnrm_ssp126_sst_monthly_onedeg_2015_2100.mat'],'sst',...
    'long_name','standard_name','units',...
    'lat','lon','time','yr');






