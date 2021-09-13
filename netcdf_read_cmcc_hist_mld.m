% Read CMIP6 netcdfs
% CMCC Hist
% MLD

clear all
close all

spath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%%
ncdisp([spath 'mlotst_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'])

%%
standard_name      = 'ocean_mixed_layer_thickness_defined_by_sigma_t';
long_name          = 'Ocean Mixed Layer Thickness Defined by Sigma T';
units              = 'm';
% online_operation   = 'average'
% cell_methods       = 'area: mean where sea time: mean'
% interval_operation = '1800 s'
% interval_write     = '1 month'
% _FillValue         = 1.000000020040877e+20
% missing_value      = 1.000000020040877e+20
% coordinates        = 'lat lon'
% description        = 'Sigma T is potential density referenced to ocean surface.'
% history            = 'none'
% cell_measures      = 'area: areacello'
%time units = 'days since 1850-01-01'

%% All years
ncid = netcdf.open([spath 'mlotst_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'],'NC_NOWRITE');
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
runs = find(yr>1965 & yr<=2015);
mld = double(mlotst(:,:,runs));
mld(mld >= 1.00e+19) = NaN;

%%
lat = latitude;
lon = longitude;

save([spath 'cmcc_hist_mld_monthly_1965_2014.mat'],'mld',...
    'long_name','standard_name','units',...
    'lat','lon','time','time_bnds','yr','runs');






