% Read CMIP6 netcdfs
% CMCC Hist
% Surface chl

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%%
ncdisp([fpath 'chlos_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'])

%%
% missing_value    = 1.000000020040877e+20
% coordinates      = 'depth lat lon'
standard_name      = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
% description      = 'sum of chlorophyll from all phytoplankton group concentrations.  In most models this is equal to chldiat+chlmisc, that is the sum of "Diatom Chlorophyll Mass Concentration" plus "Other Phytoplankton Chlorophyll Mass Concentration"'
long_name          = 'Surface Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water';
% history          = 'none'
units              = 'kg m-3';
% cell_measures    = 'area: areacello'
%Size:       362x292x1980
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'

%% All years
ncid = netcdf.open([fpath 'chlos_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_185001-201412.nc'],'NC_NOWRITE');
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
runs = find(yr>1965 & yr<=2015);
schl = double(chlos(:,:,runs));
schl(schl >= 1.00e+19) = NaN;

%%
lat = latitude;
lon = longitude;

save([fpath 'cmcc_hist_surf_chl_monthly_1965_2014.mat'],'schl',...
    'long_name','standard_name','units',...
    'lat','lon','time','time_bnds','yr','runs');






