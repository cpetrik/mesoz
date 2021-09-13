% MODIS-Aqua Chl used in GLM
% write netcdf for Jessica to calc biomes

clear all
close all

%% 
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_chl.mat']);
load([opath 'glm_obs_grid.mat'])

%% Reshape
ni=360;
nj=180;
chl = reshape(glmobschl,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);
lat = lat_g(1,:);
lon = lon_g(:,1);
time = 1:12;
nt = length(time);

%% nans to num
chl(isnan(chl(:))) = 1e20;

%% chl
file_ischl = [opath 'glm_obs_surf_chl_modisaqua.nc'];

ncidZM = netcdf.create(file_ischl,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'month',nt);

vidlat = netcdf.defVar(ncidZM,'Lat','double',lat_dim);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'Lon','double',lon_dim);
netcdf.putAtt(ncidZM,vidlon,'long_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'standard_name','lon');
netcdf.putAtt(ncidZM,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidZM,vidlon,'axis','X');

vidtZM = netcdf.defVar(ncidZM,'month','double',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','month');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','month');
netcdf.putAtt(ncidZM,vidtZM,'units','month');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');

vidbioZM = netcdf.defVar(ncidZM,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','surface chlorophyll');
netcdf.putAtt(ncidZM,vidbioZM,'units','mg/m3' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,chl);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
