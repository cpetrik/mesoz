% Read CMIP6 netcdfs
% UK-ESM Hist
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% mld
ncdisp([fpath 'mlotst_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc'])

%%
standard_name = 'ocean_mixed_layer_thickness_defined_by_sigma_t';
long_name     = 'Ocean Mixed Layer Thickness Defined by Sigma T';
% comment       = 'Sigma T is potential density referenced to ocean surface.'
units         = 'm';
% original_name = 'mo: (variable_name: mlotst)'
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'mlotst_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Get subset
% Time
yr = ((time+1)/360)+1850;
runs = find(yr>1965 & yr<=2015);
%runs = 181:780;
mld = mlotst(:,:,runs);

%%
save([spath 'ukesm_isimip_hist_mld_monthly_1965_2014.mat'],'mld',...
    'yr','runs','time','long_name','standard_name','latitude','longitude',...
    'units');





