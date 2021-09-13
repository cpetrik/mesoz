% Read CMIP6 netcdfs
% IPSL SSP 585
% Surface temp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

%%
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_tos_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_tos_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
tos(tos>=1e20) = nan;

%% Time
yr = ((time+1)/12)+1601;
tos(tos>=1e19) = NaN;
sst = tos;

%%
save([fpath 'ipsl_ssp585_sst_monthly_2015_2100.mat'],'sst','yr',...
    'long_name','standard_name','units','lat','lon',...
    'time');





