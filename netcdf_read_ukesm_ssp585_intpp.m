% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the DKRZ ISIMIP server
% Intpp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

%% intpp zall
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_ssp585_intpp_onedeg_global_monthly_2015_2100.nc'])

%%
standard_name  = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
long_name      = 'Primary Organic Carbon Production by All Types of Phytoplankton';
units          = 'mol m-2 s-1';
%missing_value = 1.000000020040877e+20;
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_ssp585_intpp_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Get subset
% Time
yr = ((time)/365)+1850;
npp = double(intpp);

%%
save([fpath 'ukesm_isimip_ssp585_npp_monthly_2015_2100.mat'],'npp',...
    'yr','time','long_name','standard_name','lat','lon','units');





