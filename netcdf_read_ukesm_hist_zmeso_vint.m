% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Meso Zoop zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-2';
missing_value = 1.000000020040877e+20;
comment       = 'vertically integrated over all ocean levels by ISIMIP data management team'
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 330
%lon: 360
%lev: 75
%zmeso: 360x180x1980
%time: 1980
%NaNs = 1.000000020040877e+20

%% Get subsets of zmeso so smaller memory
% Time
yr1 = ((time+1)/12)+1601;
runs1 = find(yr1>1950 & yr1<=2015);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs1(1)-1],[360 180 length(runs1)]);
    
end
netcdf.close(ncid);

zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;
yr = yr1;
runs = runs1;

%%
save([fpath 'ukesm_isimip_hist_zmeso_vint_monthly_1950_2014.mat'],'zmeso_vint',...
    'yr','runs','long_name','standard_name','units','lat','lon');






