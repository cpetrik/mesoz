% CMIP6 output
% 200m integrations
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('Hist_ts_chl_clim_biomes_areaw.mat')
load('Hist_ts_sst_clim_biomes_areaw.mat')
load('Hist_ts_mesoz_clim_biomes_areaw.mat')
load('glm100_output_by_biom_v2.mat')

%% Plots

cb=[34/255 136/255 51/255;...   %green
    %238/255 119/255 51/255;...  %orange
    %0/255 153/255 136/255;...   %teal
    153/255 153/255 51/255;...   %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];                     %black

set(groot,'defaultAxesColorOrder',cb);

%% Std devs as shadded region
%create continuous x value array for plotting
mo = 1:12;

X=[mo fliplr(mo)];
%create y values for out and then back
%+/- 2 stdev
Sz=[OGmzmean+2*OGmzsd; flipud(OGmzmean-2*OGmzsd)];

%% Try using tiles
close all

f1 = figure('Units','inches','Position',[1 3 12 3]);

tiledlayout(1,3, 'TileSpacing', 'compact')
% SST --------------------------
nexttile %subplot(4,3,3)
plot(mo,Ct,'LineWidth',1.5); hold on;
plot(mo,Mt,'LineWidth',1.5); hold on;
plot(mo,Nt,'LineWidth',1.5); hold on;
plot(mo,Gt,'LineWidth',1.5); hold on;
plot(mo,It,'LineWidth',1.5); hold on;
plot(mo,Ut,'LineWidth',1.5); hold on;
plot(mo,Ot,'LineWidth',1.5); hold on;
title('SST (^oC)')
%ylabel('Global')
xlim([1 12])
%text(1,16.5,'c','FontWeight','Bold','FontSize',14)

% CHL --------------------------------
nexttile %subplot(4,3,2)
plot(mo,Cc,'LineWidth',1.5); hold on;
plot(mo,Mc,'LineWidth',1.5); hold on;
plot(mo,Nc,'LineWidth',1.5); hold on;
plot(mo,Gc,'LineWidth',1.5); hold on;
plot(mo,Ic,'LineWidth',1.5); hold on;
plot(mo,Uc,'LineWidth',1.5); hold on;
plot(mo,Oc,'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
%ylabel('Global')
xlabel('NH month')
xlim([1 12])
%text(1,0.525,'b','FontWeight','Bold','FontSize',14)

% MESOZ ------------------------------
%subplot(4,3,1)
nexttile
plot(mo,Cz,'LineWidth',1.5); hold on;
plot(mo,Mz,'LineWidth',1.5); hold on;
plot(mo,Nz,'LineWidth',1.5); hold on;
plot(mo,Gz,'LineWidth',1.5); hold on;
plot(mo,Iz,'LineWidth',1.5); hold on;
plot(mo,Uz,'LineWidth',1.5); hold on;
%plot(mo,Oz,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,1),'LineWidth',1.5); hold on;
fill(X,Sz(:,1),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('zmeso (mgC m^-^2)')
ylabel('Global')
xlim([1 12])
ylim([0 1500])
%text(1,1970,'a','FontWeight','Bold','FontSize',14)

lg = legend(nexttile(3),{'CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLMM/sat'});
lg.Location = 'eastoutside';

print('-dpng',[ppath 'Hist_clim_global_means_zmeso200_schl_sst_glmm100SDs_horiz.png'])

%% Vertical
f2 = figure('Units','inches','Position',[1 3 5 9]);

%tiledlayout('vertical') - didn't work
tiledlayout(3,1, 'TileSpacing', 'compact')

% SST --------------------------
nexttile %subplot(4,3,3)
plot(mo,Ct,'LineWidth',1.5); hold on;
plot(mo,Mt,'LineWidth',1.5); hold on;
plot(mo,Nt,'LineWidth',1.5); hold on;
plot(mo,Gt,'LineWidth',1.5); hold on;
plot(mo,It,'LineWidth',1.5); hold on;
plot(mo,Ut,'LineWidth',1.5); hold on;
plot(mo,Ot,'LineWidth',1.5); hold on;
title('SST (^oC)')
%ylabel('Global')
xlim([1 12])
%text(1,16.5,'c','FontWeight','Bold','FontSize',14)

% CHL --------------------------------
nexttile %subplot(4,3,2)
plot(mo,Cc,'LineWidth',1.5); hold on;
plot(mo,Mc,'LineWidth',1.5); hold on;
plot(mo,Nc,'LineWidth',1.5); hold on;
plot(mo,Gc,'LineWidth',1.5); hold on;
plot(mo,Ic,'LineWidth',1.5); hold on;
plot(mo,Uc,'LineWidth',1.5); hold on;
plot(mo,Oc,'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
%ylabel('Global')
xlim([1 12])
%text(1,0.525,'b','FontWeight','Bold','FontSize',14)

% MESOZ ------------------------------
%subplot(4,3,1)
nexttile
plot(mo,Cz,'LineWidth',1.5); hold on;
plot(mo,Mz,'LineWidth',1.5); hold on;
plot(mo,Nz,'LineWidth',1.5); hold on;
plot(mo,Gz,'LineWidth',1.5); hold on;
plot(mo,Iz,'LineWidth',1.5); hold on;
plot(mo,Uz,'LineWidth',1.5); hold on;
%plot(mo,Oz,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,1),'LineWidth',1.5); hold on;
fill(X,Sz(:,1),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('zmeso (mgC m^-^2)')
%ylabel('Global')
xlabel('N Hem month')
xlim([1 12])
ylim([0 1500])
%text(1,1970,'a','FontWeight','Bold','FontSize',14)

lg = legend(nexttile(2),{'CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLMM/sat'});
lg.Location = 'eastoutside';

print('-dpng',[ppath 'Hist_clim_global_means_zmeso200_schl_sst_glmm100SDs_vert.png'])

