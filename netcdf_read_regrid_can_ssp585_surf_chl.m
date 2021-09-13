% Read CMIP6 CanESM5 SSP585 netcdfs
% surface chl

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';

%% surf chl
ncdisp([fpath 'chlos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'])

%%
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Surface Mass Concentration of Total Phytoplankton Expressed as Chlorophyll in Sea Water';
% comment     = 'Sum of chlorophyll from all phytoplankton group concentrations at the sea surface.  In most models this is equal to chldiat+chlmisc, that is the sum of 'Diatom Chlorophyll Mass Concentration' plus 'Other Phytoplankton Chlorophyll Mass Concentration''
units         = 'kg m-3';
% original_name = 'NCHL'
% history       = 'sum2_deptht_l1_mltby1em6'
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x180x1032
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'chlos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'],'NC_NOWRITE');
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
schl = chlos;

%%
save([fpath 'can_ssp585_surf_chl_monthly_onedeg_2015-2100.mat'],'schl',...
    'long_name','standard_name','units',...
    'lat','lon','time','yr');
