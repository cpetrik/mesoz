% Calculate different skill metrics for each ESM
% log transform biomass
% uses SH shift

clear all
close all

%%
load('skill_scores_model_obsglm_clim_raw.mat')
rskill = skill;
clear skill
load('skill_scores_model_obsglm_clim_norm.mat')
nskill = skill;
clear skill
load('skill_scores_model_obsglm_clim_anom.mat')
askill = skill;
clear skill
load('skill_scores_model_obsglm_clim_01.mat')
zskill = skill;
clear skill
load('skill_scores_model_obsglm_clim_neg11.mat')
oskill = skill;
clear skill

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey
set(groot,'defaultAxesColorOrder',cm);

ctext = {'All','Win','Spr','Sum','Fal'};

% Bar graphs
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
subplot(5,3,1)
bar(squeeze(rskill(1,:,:))')
title('Correlation coefficient')
ylabel('Raw')
set(gca,'XTickLabel',ctext)
subplot(5,3,2)
bar(squeeze(rskill(2,:,:))')
title('Root mean square error')
set(gca,'XTickLabel',ctext)
subplot(5,3,3)
bar(squeeze(rskill(11,:,:))')
title('Bias')
set(gca,'XTickLabel',ctext)

subplot(5,3,4)
bar(squeeze(nskill(1,:,:))')
ylabel('Norm')
set(gca,'XTickLabel',ctext)
subplot(5,3,5)
bar(squeeze(nskill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(5,3,6)
bar(squeeze(nskill(11,:,:))')
set(gca,'XTickLabel',ctext)

subplot(5,3,7)
bar(squeeze(askill(1,:,:))')
ylabel('Anom')
set(gca,'XTickLabel',ctext)
subplot(5,3,8)
bar(squeeze(askill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(5,3,9)
bar(squeeze(askill(11,:,:))')
set(gca,'XTickLabel',ctext)

subplot(5,3,10)
bar(squeeze(zskill(1,:,:))')
ylabel('0 to 1')
set(gca,'XTickLabel',ctext)
subplot(5,3,11)
bar(squeeze(zskill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(5,3,12)
bar(squeeze(zskill(11,:,:))')
set(gca,'XTickLabel',ctext)

subplot(5,3,13)
bar(squeeze(oskill(1,:,:))')
ylabel('-1 to 1')
set(gca,'XTickLabel',ctext)
subplot(5,3,14)
bar(squeeze(oskill(2,:,:))')
set(gca,'XTickLabel',ctext)
subplot(5,3,15)
bar(squeeze(oskill(11,:,:))')
set(gca,'XTickLabel',ctext)
print('-dpng',[figp 'corr_rmse_bias_clims_200_obsglm_all_scales.png'])

