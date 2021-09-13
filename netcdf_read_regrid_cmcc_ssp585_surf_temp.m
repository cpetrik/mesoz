% Read CMIP6 netcdfs
% CMCC SSP 585
% Surface temp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/ssp585/';

%%
ncdisp([fpath 'sst_Omon_CMCC-ESM2_ssp585_r1i1p1f1_gn_2015_2100_onedeg.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
% cell_measures = 'area: areacello'
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'sst_Omon_CMCC-ESM2_ssp585_r1i1p1f1_gn_2015_2100_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
save([fpath 'cmcc_ssp585_sst_monthly_onedeg_2015_2100.mat'],'sst',...
    'long_name','standard_name','units',...
    'lat','lon','time');






