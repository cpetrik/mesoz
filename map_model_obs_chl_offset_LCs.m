% Linear model of ESM and obs chl
% only 30S-30N
% all in units of log10 g/m3 
% see if can use offset to define LCs - forced thru zero, so no interept to
% be offset

clear all
close all

%% LM coeffs
load('coeffs_mod_obs_log10chl_tropics_only.mat')
%from log10 to raw g/m3
cff = 10.^coeffs;

% Obs LC
obsLC = log10(0.125*1e-3);
rawLC = (0.125*1e-3) + cff;
logLC = obsLC + coeffs;

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
gunits = units;
gtime = yr(runs);

gz = fliplr(gz);
pcolor(gz);

clear schl units runs yr 

%% Chl
fpath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/';
load([fpath 'seawifs_chl_ocx_growingseason_mean.mat'])

chl_ocx = chl_ocx*1e-3;
[lat_c,lon_c] = meshgrid(lat,lon);
figure
pcolor(lon_c,lat_c,chl_ocx);
figure
pcolor(chl_ocx);

%% Quantiles
qCAN(1,:) = quantile(cz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(2,:) = quantile(nz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(3,:) = quantile(gz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(4,:) = quantile(iz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(5,:) = quantile(uz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(6,:) = quantile(chl_ocx(:),[0.05 0.25 0.5 0.75 0.95]);

Tq = array2table(qCAN,'VariableNames',{'5th','25th','50th','75th','95th'},...
    'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM','obs'});
writetable(Tq,'chl_quantiles_gm3_units.csv','WriteRowNames',true)

%%
qQ(1,:) = quantile(cz(:),[0.30:0.001:0.305]);
qQ(2,:) = quantile(nz(:),[0.30:0.001:0.305]);
qQ(3,:) = quantile(gz(:),[0.30:0.001:0.305]);
qQ(4,:) = quantile(iz(:),[0.30:0.001:0.305]);
qQ(5,:) = quantile(uz(:),[0.30:0.001:0.305]);
qQ(6,:) = quantile(chl_ocx(:),[0.30:0.001:0.305]);

TQ = array2table(qQ,'VariableNames',{'300','301','302','303','304','306'},...
    'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM','obs'});
writetable(TQ,'chl_LC_quantiles_gm3_units.csv','WriteRowNames',true)

%% Maps
%add LC contour

clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

figure(1)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(cz))
cmocean('tempo')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(1,1)),'LineColor','r')
title('CAN')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(nz))
cmocean('tempo')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(2,1)),'LineColor','r')
title('CNRM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(uz))
cmocean('tempo')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(5,1)),'LineColor','r')
title('UKESM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(iz))
cmocean('tempo')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'chl (g m^-^3)')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(4,1)),'LineColor','r')
title('IPSL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(gz))
cmocean('tempo')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(3,1)),'LineColor','r')
title('GFDL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,log10(chl_ocx))
cmocean('tempo')
caxis([-5 -2])
contourm(lat_g,lon_g,log10(chl_ocx),log10(qQ(6,1)),'LineColor','r')
title('obs growing season')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_obs_chl_log.png')

%% not log10, 30th percentil
figure(2)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(cz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(1,1)),'LineColor','r')
title('CAN')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(nz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(2,1)),'LineColor','r')
title('CNRM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(uz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(5,1)),'LineColor','r')
title('UKESM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(iz))
cmocean('tempo')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'chl (g m^-^3)')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(4,1)),'LineColor','r')
title('IPSL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(gz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(3,1)),'LineColor','r')
title('GFDL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(chl_ocx))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qQ(6,1)),'LineColor','r')
title('obs growing season')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_obs_chl.png')

%% not log10, 25th percentile
figure(3)
subplot('Position',[0.01 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(cz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(1,2)),'LineColor','r')
title('CAN')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.65 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(nz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(2,2)),'LineColor','r')
title('CNRM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(uz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(5,2)),'LineColor','r')
title('UKESM')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.34 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(iz))
cmocean('tempo')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'chl (g m^-^3)')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(4,2)),'LineColor','r')
title('IPSL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(gz))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(3,2)),'LineColor','r')
title('GFDL')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.03 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,(chl_ocx))
cmocean('tempo')
caxis([0 0.00225])
contourm(lat_g,lon_g,(chl_ocx),(qCAN(6,2)),'LineColor','r')
title('obs growing season')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_obs_chl_25th.png')







