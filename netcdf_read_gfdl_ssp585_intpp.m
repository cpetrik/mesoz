% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% Intpp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% intpp zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_intpp_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
long_name     = 'Primary Organic Carbon Production by All Types of Phytoplankton';
units         = 'mol m-2 s-1';
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell volume (volcello).'
%Size:       360x180x35x1032
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_intpp_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time)/12)+1601;
npp = double(intpp);

%%
save([fpath 'gfdl_ssp585_npp_monthly_2015_2100.mat'],'npp','yr',...
    'long_name','standard_name','units','lat','lon','time');





