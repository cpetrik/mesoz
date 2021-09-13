% Read CMIP6 CanESM5 Historic netcdfs
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

%%
ncdisp([fpath 'mlotst_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'])

%%
standard_name = 'ocean_mixed_layer_thickness_defined_by_sigma_t';
long_name     = 'Ocean Mixed Layer Thickness Defined by Sigma T';
% comment       = 'Sigma T is potential density referenced to ocean surface.'
units         = 'm';
% original_name = 'mlotst'
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% time units = 'days since 1850-01-01'
% calendar   = '365_day'

%%
ncid = netcdf.open([fpath 'mlotst_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time)/365)+1850;
runs = find(yr>1965 & yr<=2015);
mld = mlotst(:,:,runs);

%%
save([spath 'can_hist_mld_monthly_1965_2014.mat'],'mld',...
    'long_name','standard_name','units',...
    'latitude','longitude','time','time_bnds','yr','runs');
