% Read CMIP6 netcdfs
% CMCC-ESM Hist
% Intpp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%%
ncdisp([fpath 'intpp_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'])

%%
long_name          = 'Primary Organic Carbon Production by All Types of Phytoplankton';
standard_name      = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
units              = 'mol m-2 s-1';
% online_operation   = 'average'
% cell_methods       = 'area: mean where sea time: mean'
% interval_operation = '1800 s'
% interval_write     = '1 month'
% _FillValue         = 1.000000020040877e+20
% missing_value      = 1.000000020040877e+20
%Size:       362x292x1980
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% All years
ncid = netcdf.open([fpath 'intpp_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time+1)/365)+1850;
runs = find(yr>1950 & yr<=2015);
npp = double(intpp(:,:,runs));

%%
save([fpath 'cmcc_hist_npp_monthly_1951_2014.mat'],'npp',...
    'long_name','standard_name','units',...
    'latitude','longitude','time','time_bnds','yr','runs');






