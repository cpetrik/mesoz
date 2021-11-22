% obsGLMM output with 100um mesh net
% write netcdf for 4z MARBL

clear all
close all

%% 
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load('glm100_obs_mesoz.mat');
load('glm100_obs_grid.mat');

%% Reshape
ni=360;
nj=180;
glmz_mo = reshape(glmobsmesoz,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);
time = 1:12;
nt = length(time);
lat = lat_g(1,:);
lon = lon_g(:,1);

%% nans to num
glmz_mo(isnan(glmz_mo(:))) = 1e20;

%% defs
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'units_vint','long_name','standard_name');

%% file
file_ischl = [opath 'obsglmm_zmeso_vint_200m_monthly_climatology.nc'];

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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name',long_name);
netcdf.putAtt(ncidZM,vidbioZM,'standard_name',standard_name);
netcdf.putAtt(ncidZM,vidbioZM,'units',units_vint);
netcdf.defVarFill(ncidZM,vidbioZM,false,1e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,glmz_mo);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
ncdisp([opath 'obsglmm_zmeso_vint_200m_monthly_climatology.nc'])

