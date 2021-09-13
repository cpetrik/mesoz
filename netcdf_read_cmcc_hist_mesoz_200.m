% Read CMIP6 netcdfs
% CMCC Hist
% Integrate top 200 m

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%% thkcello
load([fpath 'thkcello_Ofx_CMCC-ESM2_historical_r1i1p1f1_gn.mat']);

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_201001-201412.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       362x292x50x240 (last one is 60)
%Dimensions: i,j,lev,time
%time units = 'days since 1850-01-01'

%% chl time as a guide
load([fpath 'cmcc_hist_surf_chl_monthly_1965_2014.mat'],'time','yr','runs');
ctime = time;
cyr = yr;
cruns = runs;

clear time yr runs

%%
tstart = 195001:2000:201001;
tend = 196912:2000:201412;
tend = [tend 201412];

%each file has 20 years of 12 months, except last is 5 yrs of 12 months
mos = 20*12*(length(tstart) - 1) + 5*12;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
mod_yr = nan(mos,1);
mz_200 = nan(362,292,mos);

mstart = 1:240:mos;
mend = [240:240:mos mos];

%%
for t=1:length(tstart)
    
    clear latitude lev lev_bnds longitude zmeso time time_bnds i j
    clear vertices_longitude vertices_latitude ncid
    
    %%
    ncid = netcdf.open([fpath 'zmeso_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_',...
        num2str(tstart(t)),'-',num2str(tend(t)),'.nc'],'NC_NOWRITE');
    
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    
    % Get all other vars 1st
    for n = 1:(nvars)
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
    netcdf.close(ncid);
    
    %% Time
    yr1 = ((time+1)/365)+1850;
    runs1 = 1:length(yr1);
    %200 m
    z200 = find(lev <= 200);
    
    zmeso(zmeso >= 1.00e+20) = NaN;
    [ni,nj,nz,nt] = size(zmeso);
    
    %% Integrate
    zmeso = zmeso(:,:,z200,:);
    thkcello = thkcello(:,:,z200);
    zmeso_int = squeeze(nansum( ( zmeso.*repmat(thkcello,1,1,1,nt) ) ,3));
    
    %% concatenate
    mod_time(mstart(t):mend(t)) = time;
    tbnds(:,mstart(t):mend(t)) = time_bnds;
    mod_yr(mstart(t):mend(t)) = yr1;
    mz_200(:,:,mstart(t):mend(t)) = zmeso_int;
    
end

%%
runs = find(mod_yr>1965 & mod_yr<=2015);
zmeso_200 = double(mz_200(:,:,runs));

%%
lat = latitude;
lon = longitude;

%%
units_orig = units;
units_vint = 'mol m-2';

save([fpath 'cmcc_hist_zmeso_200_monthly_1965_2014.mat'],'zmeso_200','mod_yr',...
    'long_name','standard_name','units','lat','lon',...
    'lev','units_orig','units_vint');






