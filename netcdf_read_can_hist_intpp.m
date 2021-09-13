% Read CMIP6 CanESM5 Historic netcdfs
% intpp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

%%
ncdisp([fpath 'intpp_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'])

%%
standard_name = 'net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton';
long_name     = 'Primary Organic Carbon Production by All Types of Phytoplankton';
% comment       = 'Vertically integrated total primary (organic carbon) production by phytoplankton.  This should equal the sum of intpdiat+intpphymisc, but those individual components may be unavailable in some models.'
units         = 'mol m-2 s-1';
% original_name = 'PPPHY'
% history       = 'vertical_integral_sum2'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x291x1980
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'intpp_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
yr = ((time)/365)+1850;
runs = find(yr>1951 & yr<=2015);
npp = intpp(:,:,runs);

%%
save([fpath 'can_hist_npp_monthly_1951_2014.mat'],'npp',...
    'long_name','standard_name','units',...
    'latitude','longitude','time','time_bnds','yr','runs');
