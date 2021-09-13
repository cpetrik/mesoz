% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server
% Intpp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% intpp zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_intpp_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name  = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
long_name      = 'Primary Organic Carbon Production by All Types of Phytoplankton';
units          = 'mol m-2 s-1';
%missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_intpp_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
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
yr = ((time)/12)+1601;
runs = find(yr>1951 & yr<=2015);

npp = squeeze(intpp(:,:,runs));

%%
save([fpath 'ukesm_isimip_hist_npp_monthly_1950_2014.mat'],'npp',...
    'yr','runs','time','long_name','standard_name','lat','lon',...
    'units');





