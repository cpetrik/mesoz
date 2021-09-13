% Visualize COPEPOD seasonal total zooplankton carbon data
% 13 = DJF; 14 = MAM; 15 = JJA; 16 = SON

clear all
close all

%% load
load('copepod-2012_cmass-m13-qtr.mat')

%%
ulat = unique(Latitude);
ulon = unique(Longitude);
lat = -79.875:0.25:85.875;
lon = -179.875:0.25:179.875;
[lat_c,lon_c] = meshgrid(lat,lon);

figure
pcolor(lat_c)
shading flat

figure
pcolor(lon_c)
shading flat

%% interp - data are sparse, this looks crappy
zoo_c = griddata(Longitude,Latitude,TotalCarbonMassmgCm3,lon_c,lat_c);

%% check figs
figure
pcolor(lon_c,lat_c,log10(zoo_c))
shading flat
colorbar

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_c,lon_c,log10(zoo_c))
colormap('jet')
colorbar


%% loop - slow, but works

nlat = length(lat);
nlon = length(lon);

zoo = NaN*ones(nlon,nlat);
for i=1:nlon
    for j=1:nlat
        iid = find(Longitude==lon(i));
        jid = find(Latitude==lat(j));
        id = intersect(iid,jid);
        if (~isempty(id))
            zoo(i,j) = TotalCarbonMassmgCm3(id);
        end
    end
end

%%
figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_c,lon_c,log10(zoo))
colormap('jet')
colorbar
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% interp - does 1 degree look better?
lat1 = -78:85;
lon1 = -180:180;
[lat_g,lon_g] = meshgrid(lat1,lon1);
zoo_g = griddata(Longitude,Latitude,TotalCarbonMassmgCm3,lon_g,lat_g);

%% check figs
figure
pcolor(lon_g,lat_g,log10(zoo_g))
colorbar
shading flat

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(zoo_g))
colormap('jet')
colorbar
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
