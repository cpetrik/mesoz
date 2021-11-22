% Read Stromberg netcdfs
% Annual mean and monthly

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/Stromberg_x1_all/';

%% intpp zall
ncdisp([fpath 'StrombergQTR-m01_x1.nc'])

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
sz = nan*ones(360,180,12);
for mo = 1:12
    if mo<10
        ncid = netcdf.open([fpath 'StrombergQTR-m0',num2str(mo),'_x1.nc'],'NC_NOWRITE');
    else
        ncid = netcdf.open([fpath 'StrombergQTR-m',num2str(mo),'_x1.nc'],'NC_NOWRITE');
    end

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Vertically integrate top 200 m
sz(:,:,mo) = carbon_biomass * 200;

clear carbon_biomass

end

%%
units_orig = units;
units_new = 'mg-C m^-2';

save([fpath 'StrombergQTR_clim_int200_mgCm2.mat'],'sz',...
    'long_name','units_orig','lat','lon','units_new');





