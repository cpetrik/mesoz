% CMIP6 output Historic 1951-2014
% 200m integrations of zmeso
% surface chl, SST

clear all
close all

fpath = '/Volumes/MIP/Fish-MIP/CMIP6/';

%% CAN ---------------------------------------------------------------
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso_200_monthly_1950_2014.mat'],'latitude','longitude',...
    'zmeso_200','mod_yr');

time = mod_yr;
file_czoo = [cpath 'zmeso200_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_1951-2014.nc'];
[ni,nj,nt] = size(zmeso_200);

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

netcdf.putVar(ncidZM,vidlat,latitude);
netcdf.putVar(ncidZM,vidlon,longitude);
netcdf.putVar(ncidZM,vidbioZM,zmeso_200);
netcdf.putVar(ncidZM,vidtZM,time);

netcdf.close(ncidZM);

%%
clear zmeso_200 latitude longitude time

%% CNRM ---------------------------------------------------------------
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon',...
    'zmeso_200','yr1','yr2');

file_nzoo = [npath 'zmeso200_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_1950-2014.nc'];

time = [yr1; yr2];
[ni,nj,nt] = size(zmeso_200);

%% CNRM zoo
ncidZM = netcdf.create(file_nzoo,'netcdf4');

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

%%
ncdisp(file_nzoo)




