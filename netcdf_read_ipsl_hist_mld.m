% Read CMIP6 netcdfs
% IPSL Hist
% MLD

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

%% intpp zall
ncdisp([fpath 'mlotst_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc'])

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
% time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'mlotst_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time+1)/365)+1850;
runs = find(yr>1965 & yr<=2015);
mld = mlotst(:,:,runs);

%%
save([spath 'ipsl_hist_mld_monthly_1965_2014.mat'],'mld','yr',...
    'long_name','standard_name','units','nav_lat','nav_lon',...
    'runs','time');





