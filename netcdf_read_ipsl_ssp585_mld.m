% Read CMIP6 netcdfs
% IPSL SSP 585
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

%% MLD
ncdisp([fpath 'mlotst_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-210012.nc'])

%%
standard_name      = 'ocean_mixed_layer_thickness_defined_by_sigma_t';
long_name          = 'Ocean Mixed Layer Thickness Defined by Sigma T';
units              = 'm';
% online_operation   = 'average'
% cell_methods       = 'area: mean where sea time: mean'
% interval_operation = '2700 s'
% interval_write     = '1 month'
% cell_measures      = 'area: areacello'
% _FillValue         = 1.000000020040877e+20
% missing_value      = 1.000000020040877e+20
% coordinates        = 'nav_lat nav_lon'
% description        = 'Sigma T is potential density referenced to ocean surface.'                       
% time units = 'days since 2015-01-01 00:00:00'

%%
ncid = netcdf.open([fpath 'mlotst_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-210012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time+1)/365)+2015;
runs = find(yr>2051);% & yr<=2100);
mld = mlotst(:,:,runs);

%%
save([spath 'ipsl_ssp585_mld_monthly_2051_2100.mat'],'mld','yr',...
    'long_name','standard_name','units','nav_lat','nav_lon',...
    'runs','time');





