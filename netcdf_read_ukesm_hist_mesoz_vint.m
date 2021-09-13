% Read CMIP6 netcdfs
% UK-ESM Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
original_name = 'mo: (variable_name: ZME_E3T) * (C_TO_N_RATIO_ZOO: 5.625) / (variable_name: thkcello)';
original_units= 'mmol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x330x75x180
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc'],'NC_NOWRITE');

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
%zmeso: 360x330x75x600
%time: 600
%NaNs = 1.000000020040877e+20

%% Get subsets of zmeso so smaller memory
% Time
yr1 = ((time+1)/360)+1850;
runs1 = 1:length(yr1);
z1 = 1:25;

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 330 length(z1) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
% Integrate top 100 m & subset time
zmeso1_1 = squeeze(nansum(zmeso,3));
clear zmeso

%%
z2 = 26:50;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 330 length(z2) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
% Integrate 
zmeso1_2 = squeeze(nansum(zmeso,3));
clear zmeso

zmeso11 = zmeso1_1 + zmeso1_2;
clear zmeso1_1 zmeso1_2

%%
z3 = 51:75;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 330 length(z3) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);
% Integrate 
zmeso1_3 = squeeze(nansum(zmeso,3));
clear zmeso

zmeso1 = zmeso11 + zmeso1_3;
clear zmeso11 zmeso1_3

%% Clear earlier
clear latitude lev lev_bnds longitude time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

%%
ncid = netcdf.open([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 330
%lon: 360
%lev: 75
%zmeso: 360x330x75x180
%time: 180
%NaNs = 1.000000020040877e+20

%% Time
yr2 = ((time+1)/360)+1850;

zmeso2 = squeeze(nansum(zmeso,3));

clear zmeso

%%
zmeso_vint = cat(3,zmeso1,zmeso2);

%%
save([fpath 'ukesm_hist_zmeso_vint_monthly_1950_2014.mat'],'zmeso_vint','yr1',...
    'yr2','long_name','standard_name','units','latitude','longitude',...
    'lev','lev_bnds');






