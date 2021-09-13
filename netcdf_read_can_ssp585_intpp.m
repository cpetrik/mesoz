% Read CMIP6 CanESM5 Historic netcdfs
% intpp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';

%%
ncdisp([fpath 'intpp_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012.nc'])

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
% Size:       360x291x1032
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'intpp_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012.nc'],'NC_NOWRITE');
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
npp = double(intpp);

%%
save([fpath 'can_ssp585_npp_monthly_2015_2100.mat'],'npp',...
    'long_name','standard_name','units',...
    'latitude','longitude','time','time_bnds','yr');
