% Linear regressions of zmeso biomass with surf chl
% Annual climatology, biomes and global
% Historic and SSP585
% GLMM100, Stromberg, COPEPOD obs

clear all
close all

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% Delta zmeso
load('delta_log10_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes_mgC.mat');

%biomass was on log10 scale before taking difference
log_diff = (zmeans_diff);

%% Zmeso-chl relationships
load('Hist_coeffs_mod_chl_global_obsglm100_strom_cope_mgC.mat');

%% predict mesozoo over chl values
% Models have a very large range
cmin = log10(0.01);
cmax = log10(25);
chl = logspace(cmin,cmax);

hzoo0 = Intercept + Lchl .* log10(chl);

%% color same as Taylor diagrams
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

%% linear regression of hist slope against vs. future diff
min(Lchl)
max(Lchl)

xmin = -0.35;
xmax = 3.25;
xpred = [xmin:0.05:xmax]';

mdl0 = fitlm(Lchl(1:6),log_diff(:,1));
ypred0 = predict(mdl0,xpred);
r_mdl0 = mdl0.Rsquared.Ordinary;
p_mdl0 = mdl0.Coefficients.pValue(2);

noCmdl0 = fitlm(Lchl(2:6),log_diff(2:6,1));
noCypred0 = predict(noCmdl0,xpred);
r_noCmdl0 = noCmdl0.Rsquared.Ordinary;
p_noCmdl0 = noCmdl0.Coefficients.pValue(2);

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% Figures
x = (-3:0.5:6);
xs = -3:0.5:6;

%glmm100
Oy0 = Lchl(7)*ones(size(xs));

%stromberg
Sy0 = Lchl(8)*ones(size(xs));

%copepod
CSy0 = Lchl(9)*ones(size(xs));

CMy0= Lchl(10)*ones(size(xs));

%% Histogram and constrained projection
load('ECdistr_fig_Pdiff_Raw_hist_ssp585_tsmeans_zmeso200.mat')

%% global SPECIFY SUBPLOT POSITIONS
%f1 = figure('Units','inches','Position',[1 1 9 12]);
figure(1)
% Global
subplot('Position',[0.1 0.6 0.3 0.35])
plot(log10(chl),hzoo0,'LineWidth',1.5); hold on;
ylim([-2 5])
xlim([round(cmin) round(cmax)])
ylabel({'Global','log_1_0 zmeso (mgC m^-^2)'})
xlabel('log_1_0 chl (mg m^-^3)')
%title('Historic relationship')
%lg  = legend({'CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLMM','obsSM','obsCS','obsCM'});
lg  = legend(Models);
lg.Position(1:2) = [.82 .63];
lg.AutoUpdate = 'off';
text(-2,5.4,'a','FontWeight','Bold','FontSize',14)

subplot('Position',[0.5 0.6 0.3 0.35])
plot(Oy0,xs,'k','LineWidth',1.5); hold on;
plot(Sy0,xs,'color',[0.25 0.25 0.25],'LineWidth',1.5); hold on;
plot(CSy0,xs,'color',[0.5 0.5 0.5],'LineWidth',1.5); hold on;
plot(CMy0,xs,'color',[0.75 0.75 0.75],'LineWidth',1.5); hold on;
plot(xpred,ypred0,'--k','LineWidth',1.5); hold on;
for i=1:6
    plot(Lchl(i),log_diff(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0 2 -0.15 0.1])
ylabel('\Delta (log_1_0 zmeso)')
xlabel('Slope of historic relationship')
%title('Future change')
text(0.02,0.115,'b','FontWeight','Bold','FontSize',14)
text(1.0,0.085,['r^2 = ' sprintf('%2.2f',r_mdl0)])
text(1.0,0.065,['p = ' sprintf('%2.2f',p_mdl0)])

% Percent change from 1965
subplot('Position',[0.1 0.1 0.3 0.35])
plot(Myr,mF,'-.k','LineWidth',1.5); hold on;
plot(Myr,mE,'k','LineWidth',1.5); hold on;
fill(X,FCI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
fill(X,ECI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
ylabel('% \Delta zmeso')
xlabel('Year')
text(1965,6.5,'c','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
legend('full','constrained')
legend('location','southwest')

% PDFs of % change
subplot('Position',[0.5 0.1 0.3 0.35])
plot(pc,Fpdf_norm,'-.k','LineWidth',2); hold on;
plot(pc,Epdf_norm,'k','LineWidth',2); hold on;
xlabel('% \Delta zmeso')
ylabel('Probability density')
xlim([-30 16])
text(-29.5,0.084,'d','FontWeight','Bold','FontSize',14)

print('-dpng',[figp 'Hist_lm_SSP585_delta_EC_global_log10_ms_mgC.png'])
print('-djpeg',[figp 'Hist_lm_SSP585_delta_EC_global_log10_ms_mgC.jpeg'])

%% No CAN + global
figure(2)
% Global
subplot('Position',[0.1 0.6 0.3 0.35])
for i=2:10
    plot(log10(chl),hzoo0(i,:),'color',cb(i,:),'LineWidth',1.5); hold on;
end
ylim([1 4.5])
xlim([round(cmin) round(cmax)])
ylabel({'Global','log_1_0 zmeso (mgC m^-^2)'})
xlabel('log_1_0 chl (mg m^-^3)')
%title('Historic relationship')
lg  = legend({'CMCC','CNRM','GFDL','IPSL','UK','obsGLMM','obsSM','obsMO-S','obsMO-M'});
lg.Position(1:2) = [.82 .63];
lg.AutoUpdate = 'off';
text(-2,4.7,'a','FontWeight','Bold','FontSize',14)

subplot('Position',[0.5 0.6 0.3 0.35])
plot(Oy0,xs,'k','LineWidth',1.5); hold on;
plot(Sy0,xs,'color',[0.25 0.25 0.25],'LineWidth',1.5); hold on;
plot(CSy0,xs,'color',[0.5 0.5 0.5],'LineWidth',1.5); hold on;
plot(CMy0,xs,'color',[0.75 0.75 0.75],'LineWidth',1.5); hold on;
plot(xpred,noCypred0,'--k','LineWidth',1.5); hold on;
for i=2:6
    plot(Lchl(i),log_diff(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0.25 1.25 -0.14 0.0])
ylabel('\Delta (log_1_0 zmeso)')
xlabel('Slope of historic relationship')
%title('Future change')
text(0.25,0.007,'b','FontWeight','Bold','FontSize',14)
text(0.9,-0.01,['r^2 = ' sprintf('%2.2f',r_noCmdl0)])
%text(0.9,-0.02,['p = ' sprintf('%2.2f',p_noCmdl0)])

% Percent change from 1965
subplot('Position',[0.1 0.1 0.3 0.35])
plot(Myr,mF,'-.k','LineWidth',1.5); hold on;
plot(Myr,mE,'k','LineWidth',1.5); hold on;
fill(X,FCI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
fill(X,ECI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
ylabel('% \Delta zmeso')
xlabel('Year')
text(1965,6.5,'c','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
legend('full','constrained')
legend('location','southwest')

% PDFs of % change
subplot('Position',[0.5 0.1 0.3 0.35])
plot(pc,Fpdf_norm,'-.k','LineWidth',2); hold on;
plot(pc,Epdf_norm,'k','LineWidth',2); hold on;
xlabel('% \Delta zmeso')
ylabel('Probability density')
xlim([-30 16])
text(-29.5,0.084,'d','FontWeight','Bold','FontSize',14)

print('-dpng',[figp 'Hist_lm_SSP585_delta_EC_global_cope_log10_ms_noCAN_mgC.png'])
print('-djpeg',[figp 'Hist_lm_SSP585_delta_EC_global_cope_log10_ms_noCAN_mgC.jpeg'])
