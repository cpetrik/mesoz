% Read CMIP6 CanESM5 Historic netcdfs
% Integrate top 200 m

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

load([fpath 'can_grid.mat'],'thkcello')

%% Meso Zoop zall
%ncdisp([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_195101-196012.nc'])
ncdisp([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r2i1p2f1_gn_201101-201412.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
original_name = 'ZOO2';
missing_value = 1.0e+20;
%Size:       360x291x45x120
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01'

%% chl time as a guide
load([fpath 'can_hist_surf_chl_monthly_1950_2014.mat'],'time','yr','runs');
ctime = time;
cyr = yr;
cruns = runs;

clear time yr runs

%%
tstart = 195101:1000:201101;
tend = 196012:1000:201412;
tend = [tend 201412];

%each file has 10 years of 12 months, except last is 4 yrs of 12 months
mos = 10*12*(length(tstart) - 1) + 4*12;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
mod_yr = nan(mos,1);
mz_200 = nan(360,291,mos);

mstart = 1:120:mos;
mend = [120:120:mos mos];

%%
for t=1:length(tstart)
    
clear latitude lev lev_bnds longitude zmeso time time_bnds i j 
clear vertices_longitude vertices_latitude ncid

%%
ncid = netcdf.open([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_',...
    num2str(tstart(t)),'-',num2str(tend(t)),'.nc'],'NC_NOWRITE');
    
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%Time
yr1 = ((time+1)/365)+1850;
runs1 = 1:length(yr1);
%200 m
z200 = find(lev <= 200);

zmeso(zmeso >= 1.00e+20) = NaN;
[ni,nj,nz,nt] = size(zmeso);

%% Integrate 
zmeso = zmeso(:,:,z200,:);
thkcello = thkcello(:,:,z200);
zmeso_100 = squeeze(nansum( ( zmeso.*repmat(thkcello,1,1,1,nt) ) ,3));

%concatenate
mod_time(mstart(t):mend(t)) = time;
tbnds(:,mstart(t):mend(t)) = time_bnds;
mod_yr(mstart(t):mend(t)) = yr1;
mz_200(:,:,mstart(t):mend(t)) = zmeso_100;

end

%%
zmeso_200 = mz_200;
units_orig = units;
units_vint = 'mol m-2';
save([fpath 'can_hist_zmeso_200_monthly_1950_2014.mat'],'zmeso_200',...
    'long_name','standard_name','units_orig','units_vint','lev','lev_bnds',...
    'latitude','longitude','mod_time','tbnds','mod_yr');
