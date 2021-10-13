% Read Stromberg netcdfs
% Annual mean and monthly

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/Stromberg_x1_all/';

%% intpp zall
ncdisp([fpath 'StrombergQTR-m00_x1.nc'])

%%
%carbon_biomass
% Size:       360x180x1
% Dimensions: lon,lat,time
% Datatype:   double
% Attributes:
% _FillValue    = NaN
% regrid_method = 'bilinear'
units         = 'mg-C m^-3';
long_name     = 'Stromberg Carbon Biomass';

%%
ncid = netcdf.open([fpath 'StrombergQTR-m00_x1.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Vertically integrate top 200 m
zmo_all = carbon_biomass * 200;
units_orig = units;
units_new = 'mg-C m^-2';

%%
save([fpath 'StrombergQTR-m00_int200_mgCm2.mat'],'zmo_all','carbon_biomass',...
    'long_name','units_orig','lat','lon','units_new');





