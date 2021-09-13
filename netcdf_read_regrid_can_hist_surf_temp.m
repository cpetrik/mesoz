% Read CMIP6 CanESM5 Historic netcdfs
% surface temp
% Regridded by J Luo

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

%% sst
ncdisp([fpath 'tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412_onedeg.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
% comment     = 'Temperature of upper boundary of the liquid ocean, including temperatures below sea-ice and floating ice shelves.'
units         = 'degC';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x291x1980
% Dimensions: i,j,time
% time units = 'days since 1850-01-01'
% calendar   = '365_day'


%%
ncid = netcdf.open([fpath 'tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test=tos(:,:,1900);
pcolor(test)

%%
yr = ((time)/365)+1850;
runs = find(yr>1951 & yr<=2015);
sst = tos(:,:,runs);

%%
save([fpath 'can_hist_sst_monthly_onedeg_1951_2014.mat'],'sst',...
    'long_name','standard_name','units','runs',...
    'lat','lon','time','yr');

