% CMIP6 output SSP 585 2051-2100
% surf chl

clear all
close all

fpath = '/Volumes/MIP/Fish-MIP/CMIP6/chl/';

%% IPSL ---------------------------------------------------------------
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
load([ipath 'ipsl_ssp585_surf_chl_monthly_2015_2100.mat'],'lat','lon',...
    'schl','yr','long_name','units');

file_ischl = [fpath 'schl_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr_2051_2100.nc'];

ryr = yr;
tid = find(ryr>2051);
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
clear lat lon time schl

%% UKESM ---------------------------------------------------------------
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
load([upath 'ukesm_isimip_ssp585_surf_chl_monthly_2015_2100.mat'],'lat','lon',...
    'schl','yr','long_name','units')

file_uschl = [fpath 'schl_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gr_2051_2100.nc'];

ryr = yr;
tid = find(ryr>2051);
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
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
load([gpath 'gfdl_ssp585_surf_chl_monthly_2015_2100.mat'],'lat','lon',...
    'schl','yr');

file_gschl = [fpath 'schl_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_2051_2100.nc'];

ryr = yr;
tid = find(ryr>2051);
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


