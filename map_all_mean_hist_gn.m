% CMIP6 output

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso_vint_climatol_1950_2014.mat'],'zmc_all',...
    'latitude','longitude','units_vint');
load([cpath 'can_depth.mat'],'deptho');

%
cz = double(zmc_all);
cz(cz(:)<0) = 0;
cLAT = latitude;
cLON = longitude;
cdep = deptho;
cunits = units_vint;

clear zmc_all latitude longitude deptho units_vint

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso_vint_climatol_1950_2014.mat'],'zmc_all',...
    'lat','lon','units_vint');
load([npath 'cnrm_depth.mat'],'deptho');

%
nz = double(zmc_all);
nz(nz(:)<0) = 0;
nLAT = lat;
nLON = lon;
ndep = deptho;
nunits = units_vint;

clear zmc_all lat lon deptho units_vint

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_zmeso_vint_climatol_1950_2014.mat'],'zmc_all',...
    'lat','lon','units');
%load([upath 'ukesm_isimip_depth.mat']);

%
uz = double(zmc_all);
uz(uz(:)<0) = 0;
%udep = deptho;
uunits = units;
[uLAT,uLON] = meshgrid(lat,lon);

clear zmc_all lat lon deptho units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat',...
    'LAT','LON','deptho');
load([ipath 'ipsl_hist_zmeso_vint_climatol_1950_2014.mat'],'zmc_all','units');

%
iz = double(zmc_all);
iz(iz(:)<0) = 0;
iLAT = LAT;
iLON = LON;
idep = deptho;
iunits = units;

clear zmc_all LAT LON deptho units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat',...
    'LAT','LON','deptho');
load([gpath 'gfdl_hist_zmeso_vint_climatol_1950_2014.mat'],'zmc_all','units');

%
gz = double(zmc_all);
gz(gz(:)<0) = 0;
gLAT = LAT;
gLON = LON;
gdep = deptho;
gunits = units;

clear zmc_all LAT LON deptho units

%%
clatlim=[-90 90];
clonlim=[-180 180];


%% Maps
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(cLAT,cLON,(cz))
cmocean('tempo')
caxis([0 0.4])
title('CAN molC m^-^2')
%text(2.95,1.9,'CAN','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(nLAT,nLON,(nz))
cmocean('tempo')
caxis([0 0.4])
title('CNRM molC m^-^2')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(uLAT,uLON,(uz))
cmocean('tempo')
caxis([0 0.4])
title('UKESM molC m^-^2')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(iLAT,iLON,(iz))
cmocean('tempo')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'mesoz')
caxis([0 0.4])
title('IPSL molC m^-^2')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,(gz))
cmocean('tempo')
caxis([0 0.4])
title('GFDL molC m^-^2')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%subplot('Position',[0.41 0.06 0.4 0.3])
print('-dpng','Map_all_hist_clim_molC_rawunits.png')


%%

figure(2)
subplot('Position',[0.01 0.58 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(cLAT,cLON,(cz))
cmocean('tempo')
caxis([0 1e-1])
colorbar
title('CAN molC m^-^3')
%text(2.95,1.9,'CAN','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.58 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(nLAT,nLON,(nz))
cmocean('tempo')
colorbar
caxis([0 1e-1])
title('CNRM molC m^-^3')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% Quantiles
qCAN(1,:) = quantile(cz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(2,:) = quantile(nz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(3,:) = quantile(gz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(4,:) = quantile(iz(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(5,:) = quantile(uz(:),[0.05 0.25 0.5 0.75 0.95]);

Tq = array2table(qCAN,'VariableNames',{'5th','25th','50th','75th','95th'},...
    'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM'});
writetable(Tq,'zmeso_vint_quantiles_rawunits.csv','WriteRowNames',true)



