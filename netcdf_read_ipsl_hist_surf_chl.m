% Read CMIP6 netcdfs
% IPSL Hist
% Surface chl

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

%% chl zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_chl_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Chlorophyll in sea water';
units         = 'kg m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_chl_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time)/12)+1601;
runs = find(yr>1950 & yr<=2015);
z200 = olevel(1);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    chl = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1], [360 180 length(z200) length(runs)]);
end
netcdf.close(ncid);
chl(chl >= 1.00e+20) = NaN;

%% 
schl = squeeze(chl);

%%
save([fpath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'],'schl','yr',...
    'long_name','standard_name','units','lat','lon',...
    'runs','time','olevel');





