% CMIP6 output
% UKESM mesoz_vint prepared by ISIMIP team

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Hist
load([hpath 'ukesm_isimip_hist_zmeso_vint_climatol_1950_2014.mat']);
load([hpath 'ukesm_isimip_depth.mat']);

[LAT,LON] = meshgrid(lat,lon);

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
% DON'T divide by water depth for m-2 to m-3

uzD = double(zmc_DJF) * 12.01 * 1e3;% ./deptho;
uzJ = double(zmc_JJA) * 12.01 * 1e3;% ./deptho;
uzM = double(zmc_MAM) * 12.01 * 1e3;% ./deptho;
uzS = double(zmc_SON) * 12.01 * 1e3;% ./deptho;

% uzD(uzD(:)<0) = 0;
% uzJ(uzJ(:)<0) = 0;
% uzM(uzM(:)<0) = 0;
% uzS(uzS(:)<0) = 0;


%%
clatlim=[-90 90];
clonlim=[-180 180];

%% Maps
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(uzD))
cmocean('tempo')
caxis([2 4])
title('Winter ')
text(2.95,1.9,'UKESM','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(uzM))
cmocean('tempo')
caxis([2 4])
cb = colorbar('Position',[0.3 0.475 0.4 0.03],'orientation','horizontal');
xlabel(cb,'log_1_0 mesoz (mg m^-^2)')
title('Spring')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(uzJ))
cmocean('tempo')
caxis([2 4])
title('Summer')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(uzS))
cmocean('tempo')
caxis([2 4])
title('Fall')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_UKESM_isimip_hist_clim_mgCm-2.png')



