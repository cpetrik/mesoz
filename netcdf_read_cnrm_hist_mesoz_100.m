% Read CMIP6 netcdfs
% CNRM-ESM2-1 Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
original_name = 'mo: (variable_name: ZME_E3T) * (C_TO_N_RATIO_ZOO: 5.625) / (variable_name: thkcello)';
original_units= 'mmol m-3';
missing_value = 1.000000020040877e+20;
%Size:       362x294x75x180
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01'
                         
%% Early years
ncid = netcdf.open([fpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_195001-199912.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 294
%lon: 362
%lev: 75
%zmeso: 362x294x75x180
%time: 180
%NaNs = 1.000000020040877e+20

%% Get subset of zmeso
% Time
%yr = ((time+1)/12)+1850;
yr1 = ((time+1)/365)+1850;
runs1 = find(yr1>1950 & yr1<=2015);
z100 = find(lev <= 100);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 294 length(z100) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
zmeso1_100 = squeeze(nansum(zmeso,3));

%% Clear earlier
clear latitude lev lev_bnds longitude zmeso time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

%% Later
ncid = netcdf.open([fpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');

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
%zmeso: 360x330x75x180
%time: 180
%NaNs = 1.000000020040877e+20


%% Get subset of zmeso
% Time
%yr = ((time+1)/12)+1850;
yr2 = ((time+1)/365)+1850;
runs2 = 1:length(yr2);
z100 = find(lev <= 100);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs2(1)-1],[362 294 length(z100) length(runs2)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
zmeso2_100 = squeeze(nansum(zmeso,3));

%%
zmeso_100 = cat(3,zmeso1_100,zmeso2_100);

%%
save([fpath 'cnrm_hist_zmeso100_monthly_1950_2014.mat'],'zmeso_100','yr1',...
    'yr2','runs1','runs2','long_name','standard_name','units','z100',...
    'lev','lev_bounds','lat','lon');





