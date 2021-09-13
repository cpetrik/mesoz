% Read CMIP6 CanESM5 Historic netcdfs
% surface temp
% rename some things, remove others for regridding with NCO

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';

%% tos
ncdisp([fpath 'tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'])

%%
standard_name = 'sea_surface_temperature';
long_name     = 'Sea Surface Temperature';
% comment     = 'Temperature of upper boundary of the liquid ocean, including temperatures below sea-ice and floating ice shelves.'
units         = 'degC';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% coordinates   = 'latitude longitude'
% Size:       360x291x1980
% Dimensions: i,j,time
time_units = 'days since 1850-01-01';
calendar   = '365_day';


%%
ncid = netcdf.open([fpath 'tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nj,nt] = size(tos);

file_csst = [fpath 'reduc_tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'];
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
netcdf.putAtt(ncidST,vidtST,'units',time_units );
netcdf.putAtt(ncidST,vidtST,'calendar',calendar);
netcdf.putAtt(ncidST,vidtST,'axis','T');

vidbioST = netcdf.defVar(ncidST,'tos','double',[lon_dim,lat_dim,time_dim]);
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

netcdf.putVar(ncidST,vidlat,latitude);
netcdf.putVar(ncidST,vidlon,longitude);
netcdf.putVar(ncidST,vidbioST,tos);
netcdf.putVar(ncidST,vidtST,time);

netcdf.close(ncidST);

%%
ncdisp([fpath 'reduc_tos_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_185001-201412.nc'])

