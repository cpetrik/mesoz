% Read CMIP6 netcdfs
% CNRM-ESM2-1 SSP 585
% Integrate top 200 m

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-206412.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       362x294x75x600
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01'
                         

%% Early years
ncid = netcdf.open([fpath 'zmeso_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-206412.nc'],'NC_NOWRITE');

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
time1 = time;
yr1 = ((time+1)/365)+1850;
z200 = find(lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,0],[362 294 length(z200) length(time)]);
end
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'thkcello_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-206412.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% Get all other vars 1st
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[362 294 length(z200) length(time)]);
end
thkcello(thkcello >= 1.00e+20) = NaN;
netcdf.close(ncid);
netcdf.close(tcid);

%% Make sure oriented the same
whos zmeso thkcello
ztest = squeeze(zmeso(:,:,1,350));
ttest = squeeze(thkcello(:,:,1,350));

figure
pcolor(ztest)
figure
pcolor(ttest)

%% Integrate top 200 m 
zmeso1 = squeeze(nansum((zmeso.*thkcello),3));
clear zmeso thkcello 

% Clear earlier
clear latitude lev lev_bnds longitude time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

%% Later
ncid = netcdf.open([fpath 'zmeso_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_206501-210012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Time
time2 = time;
yr2 = ((time+1)/365)+1850;
%200 m
z200 = find(lev <= 200);

n = nvars;
zmeso = netcdf.getVar(ncid,n-1, [0,0,0,0],[362 294 length(z200) length(time)]);
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%% thkcello
tcid = netcdf.open([fpath 'thkcello_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_206501-210012.nc'],'NC_NOWRITE');
[ndims,tvars,ngatts,unlimdimid] = netcdf.inq(tcid);
% Get all other vars 1st
for n = 1:(tvars-1)
    varname = netcdf.inqVar(tcid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
t = tvars;
thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[362 294 length(z200) length(time)]);
netcdf.close(ncid);
thkcello(thkcello >= 1.00e+20) = NaN;


%% Make sure oriented the same
whos zmeso thkcello
ztest = squeeze(zmeso(:,:,1,30));
ttest = squeeze(thkcello(:,:,1,30));

figure
pcolor(ztest)
figure
pcolor(ttest)

%% Integrate top 200 m
zmeso2 = squeeze(nansum((zmeso.*thkcello),3));

%% Put all times together
zmeso_200 = cat(3,zmeso1,zmeso2);

clear zmeso thkcello 

%%
units_orig = units;
units_vint = 'mol m-2';

save([fpath 'cnrm_ssp585_zmeso_200_monthly_2015_2100.mat'],'zmeso_200','yr1',...
    'yr2','long_name','standard_name','units','lat','lon',...
    'lev','lev_bounds','units_orig','units_vint','time1','time2');






