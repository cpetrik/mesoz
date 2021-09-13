% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server
% Int top 200 m

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Meso Zoop zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_zmeso_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

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
z200 = find(lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 180 length(z200) length(runs1)]);
    
end
netcdf.close(ncid);

zmeso(zmeso >= 1.00e+20) = NaN;
yr = yr1;
runs = runs1;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_thkcello_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate top 200 m
zmeso_200 = squeeze(nansum((zmeso.*thkcello),3));

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'],'zmeso_200',...
    'yr','runs','long_name','standard_name','lat','lon',...
    'units_orig','units_vint','z200','lev');

%%
load([fpath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'])
test=zmeso_200(:,:,100);
[LAT,LON]=meshgrid(lat,lon);
pcolor(LON,LAT,test)






