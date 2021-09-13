% Maps of ESM and obs chl
% obs chl in mg/m3 and model in kg/m3?

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_surf_chl_monthly_onedeg_1951_2014.mat'],'schl',...
    'units','runs','yr');

cz = nanmean(schl,3) * 1e3;
cz(cz(:)<0) = 0;
cunits = units;
ctime = yr(runs);

clear schl units runs yr

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_surf_chl_monthly_onedeg_1951_2014.mat'],...
    'schl','units','runs','yr');

nz = nanmean(schl,3);
nz(nz(:)<0) = 0;
nunits = units;
ntime = yr(runs);

clear schl units runs yr

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],...
    'schl','units','runs','yr');

schl = double(schl);
uz = nanmean(schl(:,:,13:780),3) * 1e3;
uz(uz(:)<0) = 0;
uunits = units;
utime = yr(runs);

clear schl units runs yr

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat']);

iz = nanmean(schl(:,:,13:780),3);
iz(iz(:)<0) = 0;
iz = double(iz);
iunits = units;
itime = yr(runs);
[lat_g,lon_g] = meshgrid(lat,lon);

clear schl units runs yr lat lon
 
%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],...
    'schl','units','runs','yr');

gz = nanmean(schl(:,:,13:780),3) * 1e3;
gz(gz(:)<0) = 0;
gz = double(gz);
gunits = units;
gtime = yr(runs);

gz = fliplr(gz);
pcolor(gz);

clear schl units runs yr 

%% Chl
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';
opath='/Volumes/MIP/Obs_Data/Chl/';
load([opath 'globcolour_soblend_mc_x1.mat'])

chl_ocx = chl*1e-3;
[lat_c,lon_c] = meshgrid(lat,lon);

test = squeeze(chl_ocx(:,:,6));
figure
pcolor(lon_c,lat_c,test);
figure
pcolor(test);

mchl = nanmean(chl_ocx,3);

%% Quantiles
qCAN(1,:) = quantile(cz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(2,:) = quantile(nz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(3,:) = quantile(gz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(4,:) = quantile(iz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(5,:) = quantile(uz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(6,:) = quantile(mchl(:),[0.05 0.25 0.5 0.75 0.95]);

Tq = array2table(qCAN,'VariableNames',{'5th','25th','50th','75th','95th'},...
    'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM','obs'});
writetable(Tq,'chl_quantiles_gm3_units.csv','WriteRowNames',true)

%% Maps
%add LC contour

clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

cmv = colormap(viridis(11));

%% not log10
figure(1)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(cz))
colormap(cmv)
caxis([0 2.25e-3])
contourm(lat_g,lon_g,(cz),0.254e-3,'LineColor','r')
title('CAN')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(nz))
colormap(cmv)
caxis([0 2.25e-3])
contourm(lat_g,lon_g,(nz),0.164e-3,'LineColor','r')
title('CNRM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(uz))
colormap(cmv)
caxis([0 2.25e-3])
contourm(lat_g,lon_g,(uz),0.138e-3,'LineColor','r')
title('UKESM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(iz))
colormap(cmv)
caxis([0 2.25e-3])
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'chl (g m^-^3)')
% colorbar('Ticks',[0,0.25e-3,0.5e-3,1e-3,1.5e-3,2e-3],'TickLabels',...
%     [0,0.25e-3,0.5e-3,1e-3,1.5e-3,2e-3])
contourm(lat_g,lon_g,(iz),0.121e-3,'LineColor','r')
title('IPSL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(gz))
colormap(cmv)
caxis([0 2.25e-3])
contourm(lat_g,lon_g,(gz),0.268e-3,'LineColor','r')
title('GFDL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(mchl))
colormap(cmv)
caxis([0 2.25e-3])
contourm(lat_g,lon_g,(mchl),0.125e-3,'LineColor','r')
title('obs')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_obs_chl_LC.png')

%% log10
figure(2)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(cz))
colormap(viridis)
caxis([-4.5 -2.5])
contourm(lat_g,lon_g,(cz),0.254e-3,'LineColor','r')
title('CAN')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(nz))
colormap(viridis)
caxis([-4.5 -2.5])
contourm(lat_g,lon_g,(nz),0.164e-3,'LineColor','r')
title('CNRM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(uz))
colormap(viridis)
caxis([-4.5 -2.5])
contourm(lat_g,lon_g,(uz),0.138e-3,'LineColor','r')
title('UKESM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(iz))
colormap(viridis)
caxis([-4.5 -2.5])
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'chl (g m^-^3)')
contourm(lat_g,lon_g,(iz),0.121e-3,'LineColor','r')
title('IPSL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(gz))
colormap(viridis)
caxis([-4.5 -2.5])
contourm(lat_g,lon_g,(gz),0.268e-3,'LineColor','r')
title('GFDL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(mchl))
colormap(viridis)
caxis([-4.5 -2.5])
contourm(lat_g,lon_g,(mchl),0.125e-3,'LineColor','r')
title('obs')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_obs_log10chl_LC.png')





