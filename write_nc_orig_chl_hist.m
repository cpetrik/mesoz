% CMIP6 output Historic 1965-2014
% MLD

clear all
close all

fpath = '/Volumes/MIP/Fish-MIP/CMIP6/chl/';

%% IPSL ---------------------------------------------------------------
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'],'lat','lon',...
    'schl','yr','runs','long_name');

file_ischl = [fpath 'schl_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_1965_2014.nc'];

ryr = yr(runs);
tid = find(ryr>1965 & ryr<=2015);
time = ryr(tid);
chl = schl(:,:,tid);
[ni,nj,nt] = size(chl);

%% IPSL schl
ncidZM = netcdf.create(file_ischl,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

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

vidtZM = netcdf.defVar(ncidZM,'time','double',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'units','year' );
netcdf.putAtt(ncidZM,vidtZM,'calendar','365_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');

vidbioZM = netcdf.defVar(ncidZM,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name',long_name);
netcdf.putAtt(ncidZM,vidbioZM,'units','kg/m3' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,chl);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
clear lat lon time schl

%% UKESM ---------------------------------------------------------------
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],'lat','lon',...
    'schl','yr','runs','long_name','units')

file_uschl = [fpath 'schl_Omon_UKESM1-0-LL_historical_r1i1p1f2_gr_1965_2014.nc'];

ryr = yr(runs);
tid = find(ryr>1965 & ryr<=2015);
time = ryr(tid);
chl = schl(:,:,tid);
[ni,nj,nt] = size(chl);

%% U schl
ncidZM = netcdf.create(file_uschl,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

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

vidtZM = netcdf.defVar(ncidZM,'time','double',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'units','year' );
netcdf.putAtt(ncidZM,vidtZM,'calendar','365_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');

vidbioZM = netcdf.defVar(ncidZM,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name',long_name);
netcdf.putAtt(ncidZM,vidbioZM,'units',units);
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,chl);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
ncdisp(file_ischl)

%%
ncdisp(file_uschl)

%% GFDL ---------------------------------------------------------------
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],'lat','lon',...
    'schl','yr','runs');

file_gschl = [fpath 'schl_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_1965_2014.nc'];

ryr = yr(runs);
tid = find(ryr>1965 & ryr<=2015);
time = ryr(tid);
chl = schl(:,:,tid);
[ni,nj,nt] = size(chl);

%% GFDL schl
ncidZM = netcdf.create(file_gschl,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

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

vidtZM = netcdf.defVar(ncidZM,'time','double',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'units','year' );
netcdf.putAtt(ncidZM,vidtZM,'calendar','365_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');

vidbioZM = netcdf.defVar(ncidZM,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name',long_name);
netcdf.putAtt(ncidZM,vidbioZM,'units',units);
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,chl);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);


