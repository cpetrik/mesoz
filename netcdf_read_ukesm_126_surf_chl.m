% Read CMIP6 netcdfs
% UK-ESM SSP 126
% NetCDF on the DKRZ ISIMIP server
% Surface chl

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp126/';

%% chl zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_ssp126_chl_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Phytoplankton expressed as Chlorophyll in Sea Water';
units         = 'kg m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_ssp126_chl_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset
% Time
yr = ((time+1)/12)+1601;
z200 = lev(1);

%%
n = nvars;
varname = netcdf.inqVar(ncid, n-1);
chl = netcdf.getVar(ncid,n-1, [0,0,0,0],[360 180 length(z200) length(time)]);
netcdf.close(ncid);

chl(chl >= 1.00e+20) = NaN;

%% mean over time
schl = squeeze(chl);

%%
save([fpath 'ukesm_isimip_ssp126_surf_chl_monthly_2015_2100.mat'],'schl',...
    'yr','time','long_name','standard_name','lat','lon',...
    'units','lev');





