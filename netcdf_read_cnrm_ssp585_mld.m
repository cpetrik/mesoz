% Read CMIP6 netcdfs
% CNRM-ESM2-1 SSP585
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';

%%
ncdisp([fpath 'mlotst_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012.nc'])

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
% ssp585ory            = 'none'
% cell_measures      = 'area: areacello'
%time units = 'days since 1850-01-01'

%% All years
ncid = netcdf.open([fpath 'mlotst_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012.nc'],'NC_NOWRITE');
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
%runs = find(yr>2051);% & yr<=2100);
runs = 433:1032;
mld = mlotst(:,:,runs);

%%
save([spath 'cnrm_ssp585_mld_monthly_2051_2100.mat'],'mld',...
    'long_name','standard_name','units',...
    'lat','lon','time','time_bounds','yr','runs');






