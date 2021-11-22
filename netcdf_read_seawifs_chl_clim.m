% Read Seawifs chl for use with Stromberg
% Annual mean and monthly

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/Stromberg_x1_all/';
cpath='/Volumes/MIP/Obs_Data/Chl/';

%% chl 12 mo
ncdisp([cpath 'SeaWiFS_Mission_Climatology.nc'])

%%
%chl_ocx
% Size:       360x180x12
% Dimensions: lon,lat,time
% Datatype:   double

%%
ncid = netcdf.open([cpath 'SeaWiFS_Mission_Climatology.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% 
chl_mo = chl_ocx;
chl_all = nanmean(chl_ocx,3);

%%
save([cpath 'SeaWiFS_Mission_Climatology.mat'],'chl_all','chl_mo',...
    'lat','lon');
save([fpath 'SeaWiFS_Mission_Climatology.mat'],'chl_all','chl_mo',...
    'lat','lon');





