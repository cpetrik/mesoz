% Read CMIP6 netcdfs
% CMCC Hist
% Surface temp
% Regridded by J Luo

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%%
ncdisp([fpath 'sst_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_1965_2014_onedeg.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
units         = 'degC';
% cell_measures = 'area: areacello'
%Size:       360x180x600
%Dimensions: i,j,time
%time units = 'days since 1850-01-01'


%% Last 50 years
ncid = netcdf.open([fpath 'sst_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_1965_2014_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test=sst(:,:,100);
pcolor(test)


%%
save([fpath 'cmcc_hist_sst_monthly_onedeg_1965_2014.mat'],'sst',...
    'long_name','standard_name','units',...
    'lat','lon','time');
