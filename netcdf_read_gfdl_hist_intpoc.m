% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% Intpoc

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

%% intpoc zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_intpoc_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'ocean_mass_content_of_particulate_organic_matter_expressed_as_carbon';
long_name     = 'Particulate Organic Carbon Content';
units         = 'kg m-2';
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell volume (volcello).'
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_intpoc_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time)/12)+1601;
runs = find(yr>1951 & yr<=2015);

%% mean over time
int_poc = squeeze(intpoc(:,:,runs));

%%
save([fpath 'gfdl_hist_int_poc_monthly_1951_2014.mat'],'int_poc','yr',...
    'long_name','standard_name','units','lat','lon',...
    'runs');





