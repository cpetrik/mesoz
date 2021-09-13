% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% Intpp

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

%% intpp zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_intpp_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
long_name     = 'Primary Organic Carbon Production by All Types of Phytoplankton';
units         = 'mol m-2 s-1';
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell volume (volcello).'
%Size:       360x180x35x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_intpp_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time)/12)+1601;
runs = find(yr>1951 & yr<=2015);

%% mean over time
npp = squeeze(intpp(:,:,runs));

%%
save([fpath 'gfdl_hist_npp_monthly_1951_2014.mat'],'npp','yr',...
    'long_name','standard_name','units','lat','lon',...
    'runs');





