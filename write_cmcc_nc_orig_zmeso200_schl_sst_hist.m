% CMIP6 output Historic 1965-2014
% 200m integrations of zmeso
% surface chl, SST

clear all
close all

cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';

%% CMCC ---------------------------------------------------------------
load([cpath 'cmcc_hist_mld_monthly_1965_2014.mat'],'mld');
load([cpath 'cmcc_hist_sst_monthly_1965_2014.mat'],'sst');
load([cpath 'cmcc_hist_surf_chl_monthly_1965_2014.mat'],'schl','yr','runs');
load([cpath 'cmcc_hist_zmeso_200_monthly_1965_2014.mat'],'lat','lon','mod_yr',...
    'zmeso_200');

time = yr(runs);
[ni,nj,nt] = size(zmeso_200);

%% Setup netcdf path to store to
fpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
fname2 = '_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_1965_2014.nc';

file_cchl = [fpath 'schl' fname2];
file_csst = [fpath 'sst' fname2];
file_czoo = [fpath 'mesoz200' fname2];
file_cmld = [fpath 'mld' fname2];

%% Set NaNs to 1e20
mld(isnan(mld)) = 1e20;
sst(isnan(sst)) = 1e20;
schl(isnan(schl)) = 1e20;
zmeso_200(isnan(zmeso_200)) = 1e20;

%% ========================== CAN =======================================
% C chl
ncidCH = netcdf.create(file_cchl,'netcdf4');

lon_dim = netcdf.defDim(ncidCH,'x',ni);
lat_dim = netcdf.defDim(ncidCH,'y',nj);
time_dim = netcdf.defDim(ncidCH,'time',nt);

vidlat = netcdf.defVar(ncidCH,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCH,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCH,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCH,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCH,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCH,'Lon','double',[lon_dim,lat_dim]);
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
netcdf.defVarFill(ncidCH,vidbioCH,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCH,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCH,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidCH,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidCH,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidCH,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidCH);

netcdf.putVar(ncidCH,vidlat,lat);
netcdf.putVar(ncidCH,vidlon,lon);
netcdf.putVar(ncidCH,vidbioCH,schl);
netcdf.putVar(ncidCH,vidtCH,time);

netcdf.close(ncidCH);

%%
ncdisp(file_cchl)

%% C sst
ncidST = netcdf.create(file_csst,'netcdf4');

lon_dim = netcdf.defDim(ncidST,'x',ni);
lat_dim = netcdf.defDim(ncidST,'y',nj);
time_dim = netcdf.defDim(ncidST,'time',nt);

vidlat = netcdf.defVar(ncidST,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidST,vidlat,'long_name','latitude');
netcdf.putAtt(ncidST,vidlat,'standard_name','lat');
netcdf.putAtt(ncidST,vidlat,'units','degrees_north');
netcdf.putAtt(ncidST,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidST,'Lon','double',[lon_dim,lat_dim]);
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
netcdf.defVarFill(ncidST,vidbioST,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidST,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidST,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidST,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidST,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidST,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidST);

netcdf.putVar(ncidST,vidlat,lat);
netcdf.putVar(ncidST,vidlon,lon);
netcdf.putVar(ncidST,vidbioST,sst);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%% C zoo
ncidZM = netcdf.create(file_czoo,'netcdf4');

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

vidbioZM = netcdf.defVar(ncidZM,'zmeso200','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','mesozooplankton biomass integrated over top 200 m');
netcdf.putAtt(ncidZM,vidbioZM,'units','molC m-2' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.00e20);
% netcdf.putAtt(ncidZM,varid,'contact','C. Petrik');
% netcdf.putAtt(ncidZM,varid,'institution','Texas A&M University');
% netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,lat);
netcdf.putVar(ncidZM,vidlon,lon);
netcdf.putVar(ncidZM,vidbioZM,zmeso_200);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
ncdisp(file_czoo)

%% C mld
ncidMLD = netcdf.create(file_cmld,'netcdf4');

lon_dim = netcdf.defDim(ncidMLD,'x',ni);
lat_dim = netcdf.defDim(ncidMLD,'y',nj);
time_dim = netcdf.defDim(ncidMLD,'time',nt);

vidlat = netcdf.defVar(ncidMLD,'Lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidMLD,vidlat,'long_name','latitude');
netcdf.putAtt(ncidMLD,vidlat,'standard_name','lat');
netcdf.putAtt(ncidMLD,vidlat,'units','degrees_north');
netcdf.putAtt(ncidMLD,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidMLD,'Lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidMLD,vidlon,'long_name','longitude');
netcdf.putAtt(ncidMLD,vidlon,'standard_name','lon');
netcdf.putAtt(ncidMLD,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidMLD,vidlon,'axis','X');

vidtMLD = netcdf.defVar(ncidMLD,'time','double',time_dim);
netcdf.putAtt(ncidMLD,vidtMLD,'long_name','time');
netcdf.putAtt(ncidMLD,vidtMLD,'standard_name','time');
netcdf.putAtt(ncidMLD,vidtMLD,'units','year' );
netcdf.putAtt(ncidMLD,vidtMLD,'calendar','365_day');
netcdf.putAtt(ncidMLD,vidtMLD,'axis','T');

vidbioMLD = netcdf.defVar(ncidMLD,'mld','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidMLD,vidbioMLD,'long_name','Ocean Mixed Layer Thickness Defined by Sigma T');
netcdf.putAtt(ncidMLD,vidbioMLD,'units','m' );
netcdf.defVarFill(ncidMLD,vidbioMLD,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidMLD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidMLD,varid,'_FillValue',1.00e20);

netcdf.endDef(ncidMLD);

netcdf.putVar(ncidMLD,vidlat,lat);
netcdf.putVar(ncidMLD,vidlon,lon);
netcdf.putVar(ncidMLD,vidbioMLD,mld);
netcdf.putVar(ncidMLD,vidtMLD,time);

netcdf.close(ncidMLD);

