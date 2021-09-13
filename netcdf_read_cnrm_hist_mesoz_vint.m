% Read CMIP6 netcdfs
% CNRM-ESM2-1 Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
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
%zmeso: 362x294x75x600
%time: 180
%NaNs = 1.000000020040877e+20

%% Get subsets of zmeso so smaller memory
% Time
yr1 = ((time+1)/365)+1850;
runs1 = 1:length(yr1);

z1 = 1:25;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 294 length(z1) length(runs1)]);
end
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'thkcello_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_195001-199912.nc'],'NC_NOWRITE');

[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% Get all other vars 1st
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs1(1)-1],[362 294 length(z1) length(runs1)]);
end
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate top 100 m 
zmeso1_1 = squeeze(nansum((zmeso.*thkcello),3));
clear zmeso thkcello 

%% 2nd set
z2 = 26:50;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 294 length(z2) length(runs1)]);
end
zmeso(zmeso >= 1.00e+20) = NaN;
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs1(1)-1],[362 294 length(z2) length(runs1)]);
end
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate 
zmeso1_2 = squeeze(nansum((zmeso.*thkcello),3));
clear zmeso thkcello

zmeso11 = zmeso1_1 + zmeso1_2;
clear zmeso1_1 zmeso1_2

%% 3rd set
z3 = 51:75;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 294 length(z3) length(runs1)]);
end
zmeso(zmeso >= 1.00e+20) = NaN;
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs1(1)-1],[362 294 length(z3) length(runs1)]);
end
thkcello(thkcello >= 1.00e+20) = NaN;

netcdf.close(ncid);
netcdf.close(tcid);

%% Integrate 
zmeso1_3 = squeeze(nansum((zmeso.*thkcello),3));
clear zmeso thkcello

zmeso1 = zmeso11 + zmeso1_3;
clear zmeso11 zmeso1_3

% Clear earlier
clear latitude lev lev_bnds longitude time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

%% Later
ncid = netcdf.open([fpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

ncid = netcdf.open([fpath 'thkcello_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
thkcello(thkcello >= 1.00e+20) = NaN;

% Time
yr2 = ((time+1)/365)+1850;

%% Integrate all depths
zmeso2 = squeeze(nansum((zmeso.*thkcello),3));

%% Put all times together
zmeso_vint = cat(3,zmeso1,zmeso2);

%%
save([fpath 'cnrm_hist_zmeso_vint_monthly_1950_2014.mat'],'zmeso_vint','yr1',...
    'yr2','long_name','standard_name','units','lat','lon',...
    'lev','lev_bounds');

%%
units_orig = units;
units_vint = 'mol m-2';
save([fpath 'cnrm_hist_zmeso_vint_monthly_1950_2014.mat'],...
    'units_orig','units_vint','-append');





