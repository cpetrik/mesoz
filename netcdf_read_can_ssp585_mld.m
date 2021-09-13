% Read CMIP6 CanESM5 SSP netcdfs
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';

%%
ncdisp([fpath 'mlotst_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012.nc'])

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
ncid = netcdf.open([fpath 'mlotst_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012.nc'],'NC_NOWRITE');
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
runs = find(yr>2051);% & yr<=2100);
mld = mlotst(:,:,runs);

%%
save([spath 'can_ssp585_mld_monthly_2051_2100.mat'],'mld',...
    'long_name','standard_name','units',...
    'latitude','longitude','time','time_bnds','yr','runs');
