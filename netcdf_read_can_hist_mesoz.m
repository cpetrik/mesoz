% Read CMIP6 CanESM5 Historic netcdfs
% Integrate over top 100 m 

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_195101-196012.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
original_name = 'mo: (variable_name: ZME_E3T) * (C_TO_N_RATIO_ZOO: 5.625) / (variable_name: thkcello)';
original_units= 'mmol m-3';
missing_value = 1.0e+20;
%Size:       360x291x45x120
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01'

%%
tstart = 195101:1000:201101;
tend = 196012:1000:201412;
tend = [tend 201412];

%each file has 10 years of 12 months, except last is 4 yrs of 12 months
mos = 10*12*(length(tstart) - 1) + 4*12;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
mz_100 = nan(360,291,mos);

mstart = 1:120:mos;
mend = [120:120:mos mos];

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_',...
    num2str(tstart(t)),'-',num2str(tend(t)),'.nc'],'NC_NOWRITE');
    
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end

% Get subset of zmeso
%Time
yr1 = ((time+1)/365)+1850;
runs1 = 1:length(yr1);
%Top 100m
z100 = find(lev <= 100);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 291 length(z100) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

% Integrate top 100 m & subset time
zmeso_100 = squeeze(nansum(zmeso,3));

%concatenate
mod_time(mstart(t):mend(t)) = time;
tbnds(:,mstart(t):mend(t)) = time_bnds;
mz_100(:,:,mstart(t):mend(t)) = zmeso_100;

%%
clear latitude lev lev_bnds longitude zmeso time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

end

%%
zmeso_100 = mz_100;
save([fpath 'cnrm_hist_zmeso100_monthly_1950_2014.mat'],'zmeso_100',...
    'long_name','standard_name','units','z100',...
    'lev','lev_bnds','latitude','longitude','mod_time','tbnds');
