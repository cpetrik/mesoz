% CMIP6 IPSL output

clear all
close all

%% Paths
hpath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat');

CGRD = GRD;
clear GRD

%% Hist
load([hpath 'ipsl_hist_zmeso_vint_climatol_1950_2014.mat']);

%% Units
%zoo: mol C m-2

% meso zoo: from molC m-2 to gC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
% DON'T divide by water depth for m-2 to m-3

izD = double(zmc_DJF) * 12.01 * 1e3;% ./deptho;
izJ = double(zmc_JJA) * 12.01 * 1e3;% ./deptho;
izM = double(zmc_MAM) * 12.01 * 1e3;% ./deptho;
izS = double(zmc_SON) * 12.01 * 1e3;% ./deptho;


%%
clatlim=[-90 90];
clonlim=[-180 180];

%% Maps
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(izD))
cmocean('tempo')
caxis([2 4])
title('Winter ')
text(2.95,1.9,'IPSL','HorizontalAlignment','center','FontWeight','bold')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(izM))
cmocean('tempo')
caxis([2 4])
cb = colorbar('Position',[0.3 0.475 0.4 0.03],'orientation','horizontal');
xlabel(cb,'log_1_0 mesoz (mg m^-^2)')
title('Spring')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(izJ))
cmocean('tempo')
caxis([2 4])
title('Summer')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(izS))
cmocean('tempo')
caxis([2 4])
title('Fall')
print('-dpng','Map_IPSL_hist_clim_mgCm-2.png')



