% Calculate different skill metrics for each ESM
% log transform biomass
% uses SH shift

clear all
close all

%%
load('skill_scores_model_obsglm_clim_raw.mat')
rskill = skill;
clear skill

load('skill_scores_model_obsglm_clim_neg11.mat')
oskill = skill;
clear skill

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
% cm=[0 0.7 0;...   %g
%     0 0 0.75;...  %b
%     0.5 0 1;...   %purple
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
%     0.35 0.35 0.35]; %grey
cb=[34/255 136/255 51/255;...   %green
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

ctext = {'All','Win','Spr','Sum','Fal'};

%% Bar graphs
f1 = figure('Units','inches','Position',[1 3 7 6]);
subplot(3,3,1)
bar(squeeze(rskill(1,:,:))')
%BarPlotBreak([10,40,1000], 50, 60, 'RPatch', 0.85);
title('Correlation coefficient')
ylabel('Raw')
set(gca,'XTickLabel',ctext)
subplot(3,3,2)
bar(squeeze(rskill(2,:,:))')
title('Root mean square error')
set(gca,'XTickLabel',ctext)
subplot(3,3,3)
bar(squeeze(rskill(11,:,:))')
title('Bias')
set(gca,'XTickLabel',ctext)

subplot(3,3,4)
bar(squeeze(rskill(1,:,:))')
legend(simtext)
legend('location','west')
ylim([-5 -2])
ylabel('Raw')
set(gca,'XTickLabel',[],'YTickLabel',[])
subplot(3,3,5)
bar(squeeze(rskill(2,:,:))')
set(gca,'XTickLabel',ctext)
ylim([0 2.5])
subplot(3,3,6)
bar(squeeze(rskill(11,:,:))')
set(gca,'XTickLabel',ctext)
ylim([-1.5 0])

subplot(3,3,7)
bar(squeeze(oskill(1,:,:))')
ylabel('-1 to 1')
set(gca,'XTickLabel',ctext)
subplot(3,3,8)
bar(squeeze(oskill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(3,3,9)
bar(squeeze(oskill(11,:,:))')
set(gca,'XTickLabel',ctext)
print('-dpng',[figp 'corr_rmse_bias_clims_200_obsglm_raw_neg11.png'])

%% Broken y
f1 = figure('Units','inches','Position',[1 3 6.5 4]);
subplot(2,3,1)
bar(squeeze(rskill(1,:,:))')
%BarPlotBreak([10,40,1000], 50, 60, 'RPatch', 0.85);
title('Correlation coefficient')
ylabel('Raw')
set(gca,'XTickLabel',ctext)
subplot(2,3,2)
bar(squeeze(rskill(2,:,:))')
title('Root mean square error')
set(gca,'XTickLabel',ctext)
subplot(2,3,3)
bar(squeeze(rskill(11,:,:))')
title('Bias')
set(gca,'XTickLabel',ctext)

subplot(2,3,4)
bar(squeeze(oskill(1,:,:))')
ylabel('-1 to 1')
set(gca,'XTickLabel',ctext)
subplot(2,3,5)
bar(squeeze(oskill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(2,3,6)
bar(squeeze(oskill(11,:,:))')
set(gca,'XTickLabel',ctext)
print('-dpng',[figp 'corr_rmse_bias_clims_200_obsglm_raw_neg11_broken.png'])

