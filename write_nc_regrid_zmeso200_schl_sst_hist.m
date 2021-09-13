% CMIP6 output Historic 1951-2014
% 200m integrations of zmeso
% surface chl, SST
% write netcdf files of regridded ESM output for Ryan

clear all
close all

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[Lat1,Lon1] = meshgrid(lat1,lon1);

%% load
load('/Volumes/MIP/Fish-MIP/CMIP6/Gridded_v2_hist_zmeso200_schl_sst_1951_2014.mat');

%% Setup netcdf path to store to
fpath = '/Volumes/MIP/Fish-MIP/CMIP6/';
fname2 = '_hist_monthly_regridded_1951_2014.nc';

file_cchl = [fpath 'canoe_schl' fname2];
file_csst = [fpath 'canoe_sst' fname2];
file_czoo = [fpath 'canoe_mesoz' fname2];

file_nchl = [fpath 'cnrm_schl' fname2];
file_nsst = [fpath 'cnrm_sst' fname2];
file_nzoo = [fpath 'cnrm_mesoz' fname2];

file_gchl = [fpath 'gfdl_schl' fname2];
file_gsst = [fpath 'gfdl_sst' fname2];
file_gzoo = [fpath 'gfdl_mesoz' fname2];

file_ichl = [fpath 'ipsl_schl' fname2];
file_isst = [fpath 'ipsl_sst' fname2];
file_izoo = [fpath 'ipsl_mesoz' fname2];

file_uchl = [fpath 'ukesm_schl' fname2];
file_usst = [fpath 'ukesm_sst' fname2];
file_uzoo = [fpath 'ukesm_mesoz' fname2];

[ni,nj,nt] = size(zG);

time = yr_gt;

%% ========================== CAN =======================================
% C chl
ncidCH = netcdf.create(file_cchl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'lon',ni);
lat_dim = netcdf.defDim(ncidCH,'lat',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCH,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCH,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCH,vidlon,'axis','X');

vidtCH = netcdf.defVar(ncidCH,'time','double',time_dim);
netcdf.putAtt(ncidCH,vidtCH,'long_name','time');
netcdf.putAtt(ncidCH,vidtCH,'standard_name','time');
netcdf.putAtt(ncidCH,vidtCH,'units','year' );
netcdf.putAtt(ncidCH,vidtCH,'calendar','365_day');
netcdf.putAtt(ncidCH,vidtCH,'axis','T');

vidbioCH = netcdf.defVar(ncidCH,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCH,vidbioCH,'long_name','surface chlorophyll concentration');
netcdf.putAtt(ncidCH,vidbioCH,'units','kg m-3' );
%netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,Lat1);
netcdf.putVar(ncidCH,vidlon,Lon1);
netcdf.putVar(ncidCH,vidbioCH,cC);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%%
ncdisp(file_cchl)

%% C sst
ncidST = netcdf.create(file_csst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'lon',ni);
lat_dim = netcdf.defDim(ncidST,'lat',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlon,'long_name','longitude');
netcdf.putAtt(ncidST,vidlon,'standard_name','lon');
netcdf.putAtt(ncidST,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidST,vidlon,'axis','X');

vidtST = netcdf.defVar(ncidST,'time','double',time_dim);
netcdf.putAtt(ncidST,vidtST,'long_name','time');
netcdf.putAtt(ncidST,vidtST,'standard_name','time');
netcdf.putAtt(ncidST,vidtST,'units','year' );
netcdf.putAtt(ncidST,vidtST,'calendar','365_day');
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'sst','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidST,vidbioST,'long_name','sea surface temperature');
netcdf.putAtt(ncidST,vidbioST,'units','degrees C' );
%netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,Lat1);
netcdf.putVar(ncidST,vidlon,Lon1);
netcdf.putVar(ncidST,vidbioST,tC);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% C zoo
ncidZM = netcdf.create(file_czoo,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'lon',ni);
lat_dim = netcdf.defDim(ncidZM,'lat',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','mgC m-2' );
%netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,Lat1);
netcdf.putVar(ncidZM,vidlon,Lon1);
netcdf.putVar(ncidZM,vidbioZM,zC);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%% ========================== CNRM =======================================
% chl
ncidCH = netcdf.create(file_nchl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'lon',ni);
lat_dim = netcdf.defDim(ncidCH,'lat',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCH,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCH,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCH,vidlon,'axis','X');

vidtCH = netcdf.defVar(ncidCH,'time','double',time_dim);
netcdf.putAtt(ncidCH,vidtCH,'long_name','time');
netcdf.putAtt(ncidCH,vidtCH,'standard_name','time');
netcdf.putAtt(ncidCH,vidtCH,'units','year' );
netcdf.putAtt(ncidCH,vidtCH,'calendar','365_day');
netcdf.putAtt(ncidCH,vidtCH,'axis','T');

vidbioCH = netcdf.defVar(ncidCH,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCH,vidbioCH,'long_name','surface chlorophyll concentration');
netcdf.putAtt(ncidCH,vidbioCH,'units','kg m-3' );
%netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,Lat1);
netcdf.putVar(ncidCH,vidlon,Lon1);
netcdf.putVar(ncidCH,vidbioCH,cN);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%% sst
ncidST = netcdf.create(file_nsst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'lon',ni);
lat_dim = netcdf.defDim(ncidST,'lat',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlon,'long_name','longitude');
netcdf.putAtt(ncidST,vidlon,'standard_name','lon');
netcdf.putAtt(ncidST,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidST,vidlon,'axis','X');

vidtST = netcdf.defVar(ncidST,'time','double',time_dim);
netcdf.putAtt(ncidST,vidtST,'long_name','time');
netcdf.putAtt(ncidST,vidtST,'standard_name','time');
netcdf.putAtt(ncidST,vidtST,'units','year' );
netcdf.putAtt(ncidST,vidtST,'calendar','365_day');
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'sst','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidST,vidbioST,'long_name','sea surface temperature');
netcdf.putAtt(ncidST,vidbioST,'units','degrees C' );
%netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,Lat1);
netcdf.putVar(ncidST,vidlon,Lon1);
netcdf.putVar(ncidST,vidbioST,tN);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% zoo
ncidZM = netcdf.create(file_nzoo,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'lon',ni);
lat_dim = netcdf.defDim(ncidZM,'lat',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','mgC m-2' );
%netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,Lat1);
netcdf.putVar(ncidZM,vidlon,Lon1);
netcdf.putVar(ncidZM,vidbioZM,zN);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%% ========================== GFDL =======================================
% chl
ncidCH = netcdf.create(file_gchl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'lon',ni);
lat_dim = netcdf.defDim(ncidCH,'lat',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCH,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCH,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCH,vidlon,'axis','X');

vidtCH = netcdf.defVar(ncidCH,'time','double',time_dim);
netcdf.putAtt(ncidCH,vidtCH,'long_name','time');
netcdf.putAtt(ncidCH,vidtCH,'standard_name','time');
netcdf.putAtt(ncidCH,vidtCH,'units','year' );
netcdf.putAtt(ncidCH,vidtCH,'calendar','365_day');
netcdf.putAtt(ncidCH,vidtCH,'axis','T');

vidbioCH = netcdf.defVar(ncidCH,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCH,vidbioCH,'long_name','surface chlorophyll concentration');
netcdf.putAtt(ncidCH,vidbioCH,'units','kg m-3' );
%netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,Lat1);
netcdf.putVar(ncidCH,vidlon,Lon1);
netcdf.putVar(ncidCH,vidbioCH,cG);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%% sst
ncidST = netcdf.create(file_gsst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'lon',ni);
lat_dim = netcdf.defDim(ncidST,'lat',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlon,'long_name','longitude');
netcdf.putAtt(ncidST,vidlon,'standard_name','lon');
netcdf.putAtt(ncidST,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidST,vidlon,'axis','X');

vidtST = netcdf.defVar(ncidST,'time','double',time_dim);
netcdf.putAtt(ncidST,vidtST,'long_name','time');
netcdf.putAtt(ncidST,vidtST,'standard_name','time');
netcdf.putAtt(ncidST,vidtST,'units','year' );
netcdf.putAtt(ncidST,vidtST,'calendar','365_day');
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'sst','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidST,vidbioST,'long_name','sea surface temperature');
netcdf.putAtt(ncidST,vidbioST,'units','degrees C' );
%netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,Lat1);
netcdf.putVar(ncidST,vidlon,Lon1);
netcdf.putVar(ncidST,vidbioST,tG);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% zoo
ncidZM = netcdf.create(file_gzoo,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'lon',ni);
lat_dim = netcdf.defDim(ncidZM,'lat',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','mgC m-2' );
%netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,Lat1);
netcdf.putVar(ncidZM,vidlon,Lon1);
netcdf.putVar(ncidZM,vidbioZM,zG);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%% ========================== IPSL =======================================
% chl
ncidCH = netcdf.create(file_ichl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'lon',ni);
lat_dim = netcdf.defDim(ncidCH,'lat',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCH,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCH,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCH,vidlon,'axis','X');

vidtCH = netcdf.defVar(ncidCH,'time','double',time_dim);
netcdf.putAtt(ncidCH,vidtCH,'long_name','time');
netcdf.putAtt(ncidCH,vidtCH,'standard_name','time');
netcdf.putAtt(ncidCH,vidtCH,'units','year' );
netcdf.putAtt(ncidCH,vidtCH,'calendar','365_day');
netcdf.putAtt(ncidCH,vidtCH,'axis','T');

vidbioCH = netcdf.defVar(ncidCH,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCH,vidbioCH,'long_name','surface chlorophyll concentration');
netcdf.putAtt(ncidCH,vidbioCH,'units','kg m-3' );
%netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,Lat1);
netcdf.putVar(ncidCH,vidlon,Lon1);
netcdf.putVar(ncidCH,vidbioCH,cI);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%% sst
ncidST = netcdf.create(file_isst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'lon',ni);
lat_dim = netcdf.defDim(ncidST,'lat',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlon,'long_name','longitude');
netcdf.putAtt(ncidST,vidlon,'standard_name','lon');
netcdf.putAtt(ncidST,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidST,vidlon,'axis','X');

vidtST = netcdf.defVar(ncidST,'time','double',time_dim);
netcdf.putAtt(ncidST,vidtST,'long_name','time');
netcdf.putAtt(ncidST,vidtST,'standard_name','time');
netcdf.putAtt(ncidST,vidtST,'units','year' );
netcdf.putAtt(ncidST,vidtST,'calendar','365_day');
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'sst','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidST,vidbioST,'long_name','sea surface temperature');
netcdf.putAtt(ncidST,vidbioST,'units','degrees C' );
%netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,Lat1);
netcdf.putVar(ncidST,vidlon,Lon1);
netcdf.putVar(ncidST,vidbioST,tI);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% zoo
ncidZM = netcdf.create(file_izoo,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'lon',ni);
lat_dim = netcdf.defDim(ncidZM,'lat',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','mgC m-2' );
%netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,Lat1);
netcdf.putVar(ncidZM,vidlon,Lon1);
netcdf.putVar(ncidZM,vidbioZM,zI);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%% ========================== UKESM =======================================
% chl
ncidCH = netcdf.create(file_uchl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'lon',ni);
lat_dim = netcdf.defDim(ncidCH,'lat',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCH,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCH,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCH,vidlon,'axis','X');

vidtCH = netcdf.defVar(ncidCH,'time','double',time_dim);
netcdf.putAtt(ncidCH,vidtCH,'long_name','time');
netcdf.putAtt(ncidCH,vidtCH,'standard_name','time');
netcdf.putAtt(ncidCH,vidtCH,'units','year' );
netcdf.putAtt(ncidCH,vidtCH,'calendar','365_day');
netcdf.putAtt(ncidCH,vidtCH,'axis','T');

vidbioCH = netcdf.defVar(ncidCH,'schl','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCH,vidbioCH,'long_name','surface chlorophyll concentration');
netcdf.putAtt(ncidCH,vidbioCH,'units','kg m-3' );
%netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,Lat1);
netcdf.putVar(ncidCH,vidlon,Lon1);
netcdf.putVar(ncidCH,vidbioCH,cU);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%% sst
ncidST = netcdf.create(file_usst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'lon',ni);
lat_dim = netcdf.defDim(ncidST,'lat',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlon,'long_name','longitude');
netcdf.putAtt(ncidST,vidlon,'standard_name','lon');
netcdf.putAtt(ncidST,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidST,vidlon,'axis','X');

vidtST = netcdf.defVar(ncidST,'time','double',time_dim);
netcdf.putAtt(ncidST,vidtST,'long_name','time');
netcdf.putAtt(ncidST,vidtST,'standard_name','time');
netcdf.putAtt(ncidST,vidtST,'units','year' );
netcdf.putAtt(ncidST,vidtST,'calendar','365_day');
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'sst','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidST,vidbioST,'long_name','sea surface temperature');
netcdf.putAtt(ncidST,vidbioST,'units','degrees C' );
%netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,Lat1);
netcdf.putVar(ncidST,vidlon,Lon1);
netcdf.putVar(ncidST,vidbioST,tU);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% zoo
ncidZM = netcdf.create(file_uzoo,'netcdf4');

lon_dim = netcdf.defDim(ncidZM,'lon',ni);
lat_dim = netcdf.defDim(ncidZM,'lat',nj);
time_dim = netcdf.defDim(ncidZM,'time',nt);

vidlat = netcdf.defVar(ncidZM,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','lat');
netcdf.putAtt(ncidZM,vidlat,'units','degrees_north');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidZM,'lon','double',[lon_dim,lat_dim]);
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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','mgC m-2' );
%netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,Lat1);
netcdf.putVar(ncidZM,vidlon,Lon1);
netcdf.putVar(ncidZM,vidbioZM,zU);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%% Just Lat & Lon
[ni,nj] = size(Lon1);

ncidLL = netcdf.create([fpath 'Lat_Lon_onedeg_grid_vec.nc'],'netcdf4');

lon_dim = netcdf.defDim(ncidLL,'lon',ni);
lat_dim = netcdf.defDim(ncidLL,'lat',nj);

vidlat = netcdf.defVar(ncidLL,'lat','double',[lat_dim]);
netcdf.putAtt(ncidLL,vidlat,'long_name','latitude');
netcdf.putAtt(ncidLL,vidlat,'standard_name','lat');
netcdf.putAtt(ncidLL,vidlat,'units','degrees_north');
netcdf.putAtt(ncidLL,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidLL,'lon','double',[lon_dim]);
netcdf.putAtt(ncidLL,vidlon,'long_name','longitude');
netcdf.putAtt(ncidLL,vidlon,'standard_name','lon');
netcdf.putAtt(ncidLL,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidLL,vidlon,'axis','X');

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidLL,varid,'creation_date',datestr(now));
% netcdf.putAtt(ncidLL,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidLL,varid,'contact','C. Petrik');

netcdf.endDef(ncidLL);

netcdf.putVar(ncidLL,vidlat,lat1);
netcdf.putVar(ncidLL,vidlon,lon1);

netcdf.close(ncidLL);

