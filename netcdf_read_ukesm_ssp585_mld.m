% Read CMIP6 netcdfs
% UK-ESM SSP585
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

%% mld
ncdisp([fpath 'mlotst_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_205001-210012.nc'])

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
ncid = netcdf.open([fpath 'mlotst_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_205001-210012.nc'],'NC_NOWRITE');
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
runs = find(yr>2051);% & yr<=2100);
mld = mlotst(:,:,runs);

%%
save([spath 'ukesm_ssp585_mld_monthly_2051_2100.mat'],'mld',...
    'yr','runs','time','long_name','standard_name','latitude','longitude',...
    'units');





