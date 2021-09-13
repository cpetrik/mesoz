% Read in netcdf chl seasonal data from J Luo
% Latitude >= 60 N: JJA
% Latitude <= -60 S: DJF
% Latitude < 60 N & >= 30 N: JJA + MAM 
% Latitude > - 60 N & <= -30 N: SON + DJF 
% Latitude between 30 N and 30 S: all seasons
% I used seawifs mission mean ocx chlorophyll, 
% did not adjust for Southern Ocean underestimation of chlorophyll 
% (limited COPEPOD data points there)

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';

%% chl during growing season
ncdisp([fpath 'seawifs_chl_ocx_growingseason_mean.nc'])

%chl_ocx
%Size: 360x180
                         
%% 
ncid = netcdf.open([fpath 'seawifs_chl_ocx_growingseason_mean.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% save
save([fpath 'seawifs_chl_ocx_growingseason_mean.mat'],...
    'lat','lon','chl_ocx');

%%
[lat_c,lon_c] = meshgrid(lat,lon);

pcolor(lon_c,lat_c,chl_ocx)
shading flat
colorbar
caxis([0 5])

%%
nn = ~isnan(chl_ocx(:));
chl = chl_ocx(nn);
latc = lat_c(nn);
lonc = lon_c(nn);

comb(:,1) = latc;
comb(:,2) = lonc;
comb(:,3) = chl;

obsmod = array2table(comb,'VariableNames',{'Lat','Lon','chl'});

writetable(obsmod,'seawifs_chl_ocx_growingseason_mean.csv')





