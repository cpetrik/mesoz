% CMIP6 output

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_npp_monthly_1951_2014.mat'],'npp',...
    'latitude','longitude','units');

%
cz = double(npp);
cz(cz(:)<0) = 0;
cnpp = nanmean(cz,3);
cnpp(cnpp(:)>=1e20) = NaN;
cLAT = latitude;
cLON = longitude;
cunits = units;

clear npp latitude longitude units cz

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_npp_monthly_1951_2014.mat'],'npp',...
    'lat','lon','units');

%
nz = double(npp);
nz(nz(:)<0) = 0;
nnpp = nanmean(nz,3);
nnpp(nnpp(:)>=1e20) = NaN;
nLAT = lat;
nLON = lon;
nunits = units;

clear npp lat lon units nz

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_npp_monthly_1951_2014.mat'],'npp',...
    'lat','lon','units');

%
uz = double(npp);
uz(uz(:)<0) = 0;
unpp = nanmean(uz,3);
unpp(unpp(:)>=1e20) = NaN;
uunits = units;
[uLAT,uLON] = meshgrid(lat,lon);

clear npp lat lon units uz

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat',...
    'LAT','LON');
load([ipath 'ipsl_hist_npp_monthly_1951_2014.mat'],'npp','units');

%
iz = double(npp);
iz(iz(:)<0) = 0;
inpp = nanmean(iz,3);
inpp(inpp(:)>=1e20) = NaN;
iLAT = LAT;
iLON = LON;
iunits = units;

clear npp LAT LON units iz

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat',...
    'LAT','LON');
load([gpath 'gfdl_hist_npp_monthly_1951_2014.mat'],'npp','units');

%
gz = double(npp);
gz(gz(:)<0) = 0;
gnpp = nanmean(gz,3);
gnpp(gnpp(:)>=1e20) = NaN;
gLAT = LAT;
gLON = LON;
gunits = units;

clear npp LAT LON units gz

%%
clatlim=[-90 90];
clonlim=[-180 180];


%% Maps
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(cLAT,cLON,log10(cnpp))
cmocean('tempo')
caxis([-7.5 -4.5])
title('CAN molC m^-^2 s^-^1')
%text(2.95,1.9,'CAN','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(nLAT,nLON,log10(nnpp))
cmocean('tempo')
caxis([-7.5 -4.5])
title('CNRM molC m^-^2 s^-^1')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(uLAT,uLON,log10(unpp))
cmocean('tempo')
caxis([-7.5 -4.5])
title('UKESM molC m^-^2 s^-^1')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(iLAT,iLON,log10(inpp))
cmocean('tempo')
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'npp')
caxis([-7.5 -4.5])
title('IPSL molC m^-^2 s^-^1')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,log10(gnpp))
cmocean('tempo')
caxis([-7.5 -4.5])
title('GFDL molC m^-^2 s^-^1')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%subplot('Position',[0.41 0.06 0.4 0.3])
print('-dpng',[ppath 'Map_all_hist_npp_molC_rawunits_log.png'])


%% Quantiles
qCAN(1,:) = quantile(cnpp(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(2,:) = quantile(nnpp(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(3,:) = quantile(gnpp(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(4,:) = quantile(inpp(:),[0.05 0.25 0.5 0.75 0.95]);
qCAN(5,:) = quantile(unpp(:),[0.05 0.25 0.5 0.75 0.95]);

Tq = array2table(qCAN,'VariableNames',{'5th','25th','50th','75th','95th'},...
    'RowNames',{'CAN','CNRM','GFDL','IPSL','UKESM'});
writetable(Tq,'npp_quantiles_rawunits.csv','WriteRowNames',true)



