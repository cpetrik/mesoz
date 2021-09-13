% Read obs SST
% 1971-2000 from OISST
% downloaded 2/16/21

clear all
close all

fpath='/Volumes/MIP/Obs_Data/OISST/';

%% sst
ncdisp([fpath 'sst.day.mean.ltm.1971-2000.nc'])

%%
% dimensions:
% lat = 720 ;
% lon = 1440 ;
% time = UNLIMITED ; // (365 currently)
% nbnds = 2 ;

% lat:actual_range = -89.875f, 89.875f ;
% lon:actual_range = 0.125f, 359.875f ;
time_units = "days since 1800-01-01 00:00:0.0" ;
%
% sst(time, lat, lon) ;
sst_units = "degC" ;
sst_missing_value = -9.96921e+36;
sst_dataset = "NOAA High-resolution Blended Analysis" ;
% sst:actual_range = -1.8f, 34.76f ;
sst_long_name = "Daily Long Term Mean Sea Surface Temperature" ;


%%
ncid = netcdf.open([fpath 'sst.day.mean.ltm.1971-2000.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == -9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

sst(sst<-9e36) = NaN;
sst = double(sst);

%%
test = squeeze(sst(:,:,6));
pcolor(test)
shading flat
colorbar

%%
save([fpath 'sst_day_ltmean_1971_2000.mat']);

%% regrid to one degree
% shift sst lon by 180
lons = lon;
lons(lons>180) = lons(lons>180)-360;

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[Lat1,Lon1] = meshgrid(lat1,lon1);
[slat,slon] = meshgrid(double(lat),double(lons));
[olat,olon] = meshgrid(double(lat),double(lon));

%% check orientations
pcolor(Lat1); shading flat; title('Lat1');

figure
pcolor(Lon1); shading flat; title('Lon1');

figure
pcolor(slat); shading flat; title('sLat');

figure
pcolor(slon); shading flat; title('sLon');

figure
pcolor(olon); shading flat; title('oLon');

figure
pcolor(test); shading flat; title('sst');

figure
pcolor(olon,olat,test); shading flat;

%%
clatlim=[-90 90];
clonlim=[-180 180];
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(slat,slon,test)

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(olat,olon,test)

%%
% [ii,ij,it] = size(sst);
% [x,y] = size(Lat1);
% tI = NaN*ones(x,y,it);
% tG = NaN*ones(x,y,it);
% for t = 1:it
%     tI(:,:,t) = interp2(slat,slon,sst(:,:,t),Lat1,Lon1);
%     tG(:,:,t) = griddata(slat,slon,sst(:,:,t),Lat1,Lon1);
% end
% 
% %%
% testI = squeeze(tI(:,:,150));
% figure
% pcolor(testI)
% shading flat
% colorbar
% 
% testG = squeeze(tG(:,:,150));
% figure
% pcolor(testG)
% shading flat
% colorbar

%% take annual mean, then regrid
msst = nanmean(sst,3);
[x,y] = size(Lat1);
%sstI = interp2(slat,slon,msst,Lat1,Lon1);
sstG = griddata(slat,slon,msst,Lat1,Lon1);

%%
% figure
% pcolor(sstI)
% shading flat
% colorbar

figure
pcolor(sstG)
shading flat
colorbar

% figure
% pcolor(sstG-sstI)
% shading flat
% colorbar

%%
save([fpath 'sst_annual_mean_onedeg_1971_2000.mat'],'Lat1','Lon1',...
    'sstG');


