% CMIP6 COBALT output

clear all
close all

%% Paths

hpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat');

CGRD = GRD;
clear GRD

%% Hist
load([hpath 'gfdl_hist_zmeso_vint_climatol_1950_2014.mat']);

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-3 to g(WW) m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
% divide by water depth for m-2 to m-3?

gzD = double(zmc_DJF) * 12.01 * 1e3 ./deptho;
gzJ = double(zmc_JJA) * 12.01 * 1e3 ./deptho;
gzM = double(zmc_MAM) * 12.01 * 1e3 ./deptho;
gzS = double(zmc_SON) * 12.01 * 1e3 ./deptho;


%%
clatlim=[-90 90];
clonlim=[-180 180];

%% GFDL
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(gzD))
cmocean('tempo')
caxis([-1 2])
title('Winter ')
text(2.9,1.9,'GFDL','HorizontalAlignment','center','FontWeight','bold')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(gzM))
cmocean('tempo')
caxis([-1 2])
cb = colorbar('Position',[0.3 0.475 0.4 0.03],'orientation','horizontal');
xlabel(cb,'log_1_0 mesoz (mg m^-^3)')
title('Spring')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(gzJ))
cmocean('tempo')
caxis([-1 2])
title('Summer')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(gzS))
cmocean('tempo')
caxis([-1 2])
title('Fall')
print('-dpng','Map_GFDL_hist_clim_mgC.png')



