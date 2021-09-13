% Read CMIP6 netcdfs
% CNRM-ESM2-1 SSP 585
% Mesozoo in top 200 m

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';

%%
ncdisp([fpath 'zmeso200_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_2015_2100_onedeg.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-2';
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'zmeso200_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_2015_2100_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time)/365)+1850;

save([fpath 'cnrm_ssp585_zmeso200_monthly_onedeg_2015_2100.mat'],'zmeso200',...
    'long_name','standard_name','units',...
    'lat','lon','time','yr');






