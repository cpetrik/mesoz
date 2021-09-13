% Visualize COPEPOD seasonal total zooplankton carbon data
% 13 = DJF; 14 = MAM; 15 = JJA; 16 = SON

clear all
close all

%% load
load('copepod-2012_cmass-m16-qtr.mat')

%% 1/2 degree
lat = -79.5:0.5:89.5;
lon = -179.5:0.5:179.5;
[lat_c,lon_c] = meshgrid(lat,lon);

%% loop - slow, but works

nlat = length(lat);
nlon = length(lon);

zoo_c = NaN*ones(nlon,nlat);
for i=1:nlon-1
    for j=1:nlat-1
        iid = find(Longitude>=lon(i) & Longitude<lon(i+1));
        jid = find(Latitude>=lat(j) & Latitude<lat(j+1));
        id = intersect(iid,jid);
        if (~isempty(id))
            zoo_c(i,j) = mean(TotalCarbonMassmgCm3(id));
        end
    end
end

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=30;
plotmaxlon=390;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% figure
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[30 390],'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(lat_c,lon_c,log10(zoo_c))
% colormap('jet')
% colorbar
% caxis([-1 2])
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% title('Fall 1/2^o')
% print('-dpng','Total_zoo_carbon_mass_fall_half.png')

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

%% loop - slow, but works

nlat1 = length(lat1);
nlon1 = length(lon1);

zoo_g = NaN*ones(nlon1,nlat1);
for i=1:nlon1-1
    for j=1:nlat1-1
        iid = find(Longitude>=lon1(i) & Longitude<lon1(i+1));
        jid = find(Latitude>=lat1(j) & Latitude<lat1(j+1));
        id = intersect(iid,jid);
        if (~isempty(id))
            zoo_g(i,j) = mean(TotalCarbonMassmgCm3(id));
        end
    end
end

%%
figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[30 390],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(zoo_g))
colormap('jet')
colorbar
caxis([-1 2])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Fall 1^o')
print('-dpng','Total_zoo_carbon_mass_fall_one.png')

%% save
fileid = COPEPOD2012fileidentifier(1);
units = 'Total Carbon Mass mgC m-3';
save('copepod-2012_cmass_SON_gridded.mat','lat_c','lon_c','zoo_c','lat_g','lon_g',...
    'zoo_g','fileid','units');






