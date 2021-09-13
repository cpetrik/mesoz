% Read CMIP6 netcdfs
% GFDL-ESM4 SSP 585
% Surface chl

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% chl zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_chl_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water';
units         = 'kg m-3';
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello volume: volcello'
% interp_method = 'conserve_order1'
% original_name = 'chl'
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell volume (volcello).'
%Size:       360x180x35x1032
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_chl_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;
z200 = lev(1);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    chl = netcdf.getVar(ncid,n-1, [0,0,0,0],[360 180 length(z200) length(time)]);
    
end
netcdf.close(ncid);
chl(chl >= 1.00e+20) = NaN;

%% mean over time
schl = squeeze(chl);

%%
save([fpath 'gfdl_ssp585_surf_chl_monthly_2015_2100.mat'],'schl','yr',...
    'long_name','standard_name','units','lat','lon',...
    'lev');





