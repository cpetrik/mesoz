% CMIP6 output

clear all
close all

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';

%% Hist
load([hpath 'cnrm_hist_zmeso_vint_climatol_1950_2014.mat']);
load([hpath 'cnrm_depth.mat'],'deptho');

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
% DON'T divide by depth from m-2 to m-3

czD = double(zmc_DJF) * 12.01 * 1e3;% ./ deptho;
czJ = double(zmc_JJA) * 12.01 * 1e3;% ./ deptho;
czM = double(zmc_MAM) * 12.01 * 1e3;% ./ deptho;
czS = double(zmc_SON) * 12.01 * 1e3;% ./ deptho;


%%
clatlim=[-90 90];
clonlim=[-180 180];

LAT = lat;
LON = lon;

%% Maps
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(czD))
cmocean('tempo')
caxis([2 4])
title('Winter ')
text(2.95,1.9,'CNRM','HorizontalAlignment','center','FontWeight','bold')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(czM))
cmocean('tempo')
caxis([2 4])
cb = colorbar('Position',[0.3 0.475 0.4 0.03],'orientation','horizontal');
xlabel(cb,'log_1_0 mesoz (mg m^-^2)')
title('Spring')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(czJ))
cmocean('tempo')
caxis([2 4])
title('Summer')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(czS))
cmocean('tempo')
caxis([2 4])
title('Fall')
print('-dpng','Map_CNRM_hist_clim_mgCm-2.png')



