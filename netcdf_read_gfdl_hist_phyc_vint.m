% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% MLD on regular grid

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

%% phyc zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_phyc_vint_onedeg_global_monthly_1850_2014.nc'])

%%
% Size:       360x180x1980
standard_name = 'mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Phytoplankton Carbon Content';
units         = 'mol m-2';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
% interp_method = 'conserve_order1'
% original_name = 'mlotst'
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'
% time units    = 'days since 1850-01-01 00:00:00'
% calendar      = 'noleap'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_phyc_vint_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

for i = (nvars)
    varname = netcdf.inqVar(ncid, i-1);
    phycvint = netcdf.getVar(ncid,i-1);
    %eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time)/12)+1601;
runs = find(yr>1951 & yr<=2015);

%% mean over time
phyc_vint = squeeze(phycvint(:,:,runs));

%%
save([fpath 'gfdl_hist_phyc_vint_monthly_1951_2014.mat'],'phyc_vint','yr',...
    'long_name','standard_name','units','lat','lon',...
    'runs');





