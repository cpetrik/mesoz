% CMIP6 output 
% 200m integrations
%Bar plot of mean biomasses

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
ddir = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%%
load([ddir 'bar_zmeso_global_biome_means.mat']);

%% colors ?
cb=[34/255 136/255 51/255;...   %green
    153/255 153/255 51/255;...  %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0;...                   %black
    0.25 0.25 0.25;...             %dk grey
    0.50 0.50 0.50;...             % grey
    0.75 0.75 0.75];               %lt grey

set(groot,'defaultAxesColorOrder',cb);

%% Update legend
Model{9} = 'obsMO-M';
Model{10} = 'obsMO-S';

%% raw
figure(1)
bar(zmeans')
legend(Model)
legend('location','northwest')
set(gca,'XTickLabels',region)
ylabel('mean zmeso biomass (mgC m^-^2)')
print('-dpng',[ppath 'Bar_Raw_hist_means_zmeso200_revis.png'])

%% lg10
figure(2)
bar(log10(zmeans)')
legend(Model)
legend('location','eastoutside')
set(gca,'XTickLabels',region)
ylabel('mean zmeso biomass (log_1_0 mgC m^-^2)')
print('-dpng',[ppath 'Bar_log10_hist_means_zmeso200_revis.png'])

%% gC
figure(3)
bar(1e-3*zmeans')
legend(Model)
legend('location','northwest')
set(gca,'XTickLabels',region)
ylabel('mean zmeso biomass (gC m^-^2)')
print('-dpng',[ppath 'Bar_Raw_hist_means_zmeso200_revis_gC.png'])


