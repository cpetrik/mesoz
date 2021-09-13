% Read CMIP6 netcdfs
% Area on regular grid

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/';

%% 
ncdisp([fpath 'ISIMIP3b_cellarea_onedeg.nc'])

%%
standard_name = 'area';
long_name     = 'area of grid cell';
units         = 'm2';

%%
ncid = netcdf.open([fpath 'ISIMIP3b_cellarea_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% grid
[lat_a,lon_a] = meshgrid(lat,lon);

figure
pcolor(cell_area) 
shading flat
title('area')

figure
pcolor(lat_a)
shading flat
title('lat')

figure
pcolor(lon_a)
shading flat
title('lon')

%% is symmetric around eq
sum(cell_area(:) == fliplr(cell_area(:)))

%%
save([fpath 'ISIMIP3b_cellarea_onedeg.mat'],'cell_area',...
    'long_name','standard_name','units','lat','lon',...
    'lat_a','lon_a');





