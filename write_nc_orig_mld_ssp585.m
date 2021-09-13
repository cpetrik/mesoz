% CMIP6 output SSP 585 2051-2100
% MLD

clear all
close all

fpath = '/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';

%% CAN ---------------------------------------------------------------
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';
load([cpath 'can_ssp585_mld_monthly_2051_2100.mat'],'latitude','longitude',...
    'mld','yr','runs');%,'time');

%mod_yr = time(runs);
mod_yr = yr(runs);
file_cmld = [fpath 'mld_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_2051_2100.nc'];
[ni,nj,nt] = size(mld);

%% C mld
ncidZM = netcdf.create(file_cmld,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'Lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidZM,vidbioZM,'units','m' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,latitude);
netcdf.putVar(ncidZM,vidlon,longitude);
netcdf.putVar(ncidZM,vidbioZM,mld);
netcdf.putVar(ncidZM,vidtZM,mod_yr);

netcdf.close(ncidZM);

%%
clear mld latitude longitude time mod_yr

%% CNRM ---------------------------------------------------------------
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';
load([npath 'cnrm_ssp585_mld_monthly_2051_2100.mat'],'lat','lon',...
    'mld','yr','runs');

file_nmld = [fpath 'mld_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_2051_2100.nc'];

time = yr(runs);
[ni,nj,nt] = size(mld);

%% CNRM mld
ncidZM = netcdf.create(file_nmld,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'Lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidZM,vidbioZM,'units','m' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,mld);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
clear lat lon time mld

%% IPSL ---------------------------------------------------------------
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
load([ipath 'ipsl_ssp585_mld_monthly_2051_2100.mat'],'nav_lat','nav_lon',...
    'mld','yr','runs');

file_imld = [fpath 'mld_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_2051_2100.nc'];

time = yr(runs);
[ni,nj,nt] = size(mld);

%% IPSL mld
ncidZM = netcdf.create(file_imld,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'Lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidZM,vidbioZM,'units','m' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,nav_lat);
netcdf.putVar(ncidZM,vidlon,nav_lon);
netcdf.putVar(ncidZM,vidbioZM,mld);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
clear nav_lat nav_lon time mld

%% UKESM ---------------------------------------------------------------
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
load([upath 'ukesm_ssp585_mld_monthly_2051_2100.mat'],'latitude','longitude',...
    'mld','yr','runs');%,'time');

time = yr(runs);
file_umld = [fpath 'mld_Omon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_2051_2100.nc'];
[ni,nj,nt] = size(mld);

%% U mld
ncidZM = netcdf.create(file_umld,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'x',ni);
lat_dim = netcdf.defDim(ncidZM,'y',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'Lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlon,'long_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'standard_name','lon');
netcdf.putAtt(ncidZM,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidZM,vidlon,'axis','X');

vidtZM = netcdf.defVar(ncidZM,'time','double',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'units','year' );
netcdf.putAtt(ncidZM,vidtZM,'calendar','360_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');

vidbioZM = netcdf.defVar(ncidZM,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidZM,vidbioZM,'units','m' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,latitude);
netcdf.putVar(ncidZM,vidlon,longitude);
netcdf.putVar(ncidZM,vidbioZM,mld);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
ncdisp(file_cmld)

%%
ncdisp(file_nmld)

%%
ncdisp(file_imld)

%%
ncdisp(file_umld)

%% GFDL ---------------------------------------------------------------
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
load([gpath 'gfdl_ssp585_mld_monthly_2051_2100.mat'],'lat','lon',...
    'mld','yr','runs');

file_gmld = [fpath 'mld_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_2051_2100.nc'];

time = yr(runs);
[ni,nj,nt] = size(mld);

%% GFDL mld
ncidZM = netcdf.create(file_gmld,'netcdf4');

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

vidbioZM = netcdf.defVar(ncidZM,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidZM,vidbioZM,'units','m' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,mld);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);



