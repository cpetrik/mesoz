% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server
% Surface chl

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% chl zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_chl_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Phytoplankton expressed as Chlorophyll in Sea Water';
units         = 'kg m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_chl_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
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
runs = find(yr>1950 & yr<=2015);
z200 = lev(1);

%%
n = nvars;
varname = netcdf.inqVar(ncid, n-1);
chl = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
netcdf.close(ncid);

chl(chl >= 1.00e+20) = NaN;

%% mean over time
schl = squeeze(chl);

%%
save([fpath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],'schl',...
    'yr','runs','time','long_name','standard_name','lat','lon',...
    'units','lev');





