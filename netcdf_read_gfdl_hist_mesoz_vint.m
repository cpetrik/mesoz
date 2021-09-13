% Read CMIP6 netcdfs
% GFDL-ESM4 Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

%% Meso Zoop vint
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-2';
missing_value = 1.000000020040877e+20;
comment       = 'vertically integrated over all ocean levels by ISIMIP data management team';
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%% 
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);

%% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1],[360 180 length(runs)]);
    
end
netcdf.close(ncid);

zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

%%
save([fpath 'gfdl_hist_zmeso_vint_monthly_1950_2014.mat'],'zmeso_vint','yr',...
    'long_name','standard_name','units','lat','lon');





