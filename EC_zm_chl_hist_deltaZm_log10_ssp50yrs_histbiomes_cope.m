% Linear regressions of zmeso biomass with surf chl 
% Annual climatology, biomes and global
% Historic and SSP585
% GLMM100, Stromberg, COPEPOD obs

clear all
close all

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% Delta zmeso
load('delta_log10_means_areaw_hist_ssp585_50yr_zmeso200_histbiomes.mat');

%biomass was on log10 scale before taking difference
log_diff = (zmeans_diff); 

%% Zmeso-chl relationships
load('Hist_SSP585_coeffs_log10_chl_global_biomes.mat');
load('Hist_CIs_biomes_obsglm100.mat');
load('Hist_std_err_biomes_obsglm100.mat');
load('Hist_CIs_std_err_biomes_stromberg.mat');
load('Hist_coeffs_chl_global_biomes_copepod.mat')

%%
Hist.slop = [HistB1(1:6),HistB2(1:6),HistB3(1:6),HistBG(1:6)];
Hist.int = [HistA1(1:6),HistA2(1:6),HistA3(1:6),HistAG(1:6)];

SSP.slop = [SSP585B1(1:6),SSP585B2(1:6),SSP585B3(1:6),SSP585BG(1:6)];
SSP.int = [SSP585A1(1:6),SSP585A2(1:6),SSP585A3(1:6),SSP585AG(1:6)];

%%
Hist.slop(7,:) = [bOG(2:4)' bOG(1)];
Hist.slop(8,:) = [bSM(2:4)' bSM(1)];
Hist.slop(9,:) = SchlB';
Hist.slop(10,:) = MchlB';

Hist.int(7,:) = [aOG(2:4)' aOG(1)];
Hist.int(8,:) = [aSM(2:4)' aSM(1)];
Hist.int(9,:) = SchlA';
Hist.int(10,:) = MchlA';

%% predict mesozoo over chl values
% Models have a very large range
cmin = log10(0.01);
cmax = log10(25);
chl = logspace(cmin,cmax);

hzoo1 = Hist.int(:,1) + Hist.slop(:,1) .* log10(chl);
szoo1 = SSP.int(:,1) + SSP.slop(:,1) .* log10(chl);
hzoo2 = Hist.int(:,2) + Hist.slop(:,2) .* log10(chl);
szoo2 = SSP.int(:,2) + SSP.slop(:,2) .* log10(chl);
hzoo3 = Hist.int(:,3) + Hist.slop(:,3) .* log10(chl);
szoo3 = SSP.int(:,3) + SSP.slop(:,3) .* log10(chl);
hzoo4 = Hist.int(:,4) + Hist.slop(:,4) .* log10(chl);
szoo4 = SSP.int(:,4) + SSP.slop(:,4) .* log10(chl);

%% color same as Taylor diagrams
cb=[34/255 136/255 51/255;...   %green
    153/255 153/255 51/255;...  %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0;...                   %black
    0.3 0.3 0.3;...             %dk grey
    0.6 0.6 0.6;...             %lt grey
    0.6 0.6 0.6];               %lt grey

set(groot,'defaultAxesColorOrder',cb);

%hmod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obs'};
hmod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','GM','SM','Cobs'};
smod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK'};

%% linear regression of hist slope against vs. future diff
min(Hist.slop(:))
max(Hist.slop(:))

xmin = -0.35;
xmax = 3.05;
xpred = [xmin:0.05:xmax]';

mdl1 = fitlm(Hist.slop(1:6,1),log_diff(:,2));
ypred1 = predict(mdl1,xpred);

mdl2 = fitlm(Hist.slop(1:6,2),log_diff(:,3));
ypred2 = predict(mdl2,xpred);

mdl3 = fitlm(Hist.slop(1:6,3),log_diff(:,4));
ypred3 = predict(mdl3,xpred);

mdl4 = fitlm(Hist.slop(1:6,4),log_diff(:,1));
ypred4 = predict(mdl4,xpred);

noCmdl1 = fitlm(Hist.slop(2:6,1),log_diff(2:6,2));
noCypred1 = predict(noCmdl1,xpred);

noCmdl2 = fitlm(Hist.slop(2:6,2),log_diff(2:6,3));
noCypred2 = predict(noCmdl2,xpred);

noCmdl3 = fitlm(Hist.slop(2:6,3),log_diff(2:6,4));
noCypred3 = predict(noCmdl3,xpred);

noCmdl4 = fitlm(Hist.slop(2:6,4),log_diff(2:6,1));
noCypred4 = predict(noCmdl4,xpred);

test=mdl1.Coefficients;
tt=table2array(test);

simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLMM','obsSM','obsCOPE'};
sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% Figures
x = (-3:0.5:6);
xs = -3:0.5:6;

%glmm100
Oy1 = bOG(2)*ones(size(xs));
Oy2 = bOG(3)*ones(size(xs));
Oy3 = bOG(4)*ones(size(xs));
Oy4 = bOG(1)*ones(size(xs));

Ol1 = bOG(2)*ones(size(xs)) - 2*bOGse(2);
Ol2 = bOG(3)*ones(size(xs)) - 2*bOGse(3);
Ol3 = bOG(4)*ones(size(xs)) - 2*bOGse(4);
Ol4 = bOG(1)*ones(size(xs)) - 2*bOGse(1);

Ou1 = bOG(2)*ones(size(xs)) + 2*bOGse(2);
Ou2 = bOG(3)*ones(size(xs)) + 2*bOGse(3);
Ou3 = bOG(4)*ones(size(xs)) + 2*bOGse(4);
Ou4 = bOG(1)*ones(size(xs)) + 2*bOGse(1);

%stromberg
Sy1 = bSM(2)*ones(size(xs));
Sy2 = bSM(3)*ones(size(xs));
Sy3 = bSM(4)*ones(size(xs));
Sy4 = bSM(1)*ones(size(xs));

Sl1 = bSM(2)*ones(size(xs)) - 2*bSMse(2);
Sl2 = bSM(3)*ones(size(xs)) - 2*bSMse(3);
Sl3 = bSM(4)*ones(size(xs)) - 2*bSMse(4);
Sl4 = bSM(1)*ones(size(xs)) - 2*bSMse(1);

Su1 = bSM(2)*ones(size(xs)) + 2*bSMse(2);
Su2 = bSM(3)*ones(size(xs)) + 2*bSMse(3);
Su3 = bSM(4)*ones(size(xs)) + 2*bSMse(4);
Su4 = bSM(1)*ones(size(xs)) + 2*bSMse(1);

%copepod
CSy1 = SchlB(1)*ones(size(xs));
CSy2 = SchlB(2)*ones(size(xs));
CSy3 = SchlB(3)*ones(size(xs));
CSy4 = SchlB(4)*ones(size(xs));

CMy1 = MchlB(1)*ones(size(xs));
CMy2 = MchlB(2)*ones(size(xs));
CMy3 = MchlB(3)*ones(size(xs));
CMy4 = MchlB(4)*ones(size(xs));

%% 3 biomes + global SPECIFY SUBPLOT POSITIONS
f1 = figure('Units','inches','Position',[1 1 9 12]);
% Global
%subplot(4,4,1)
subplot('Position',[0.1 0.765 0.25 0.2])
plot(log10(chl),hzoo4,'LineWidth',1.5); hold on;
ylim([-6 4])
xlim([round(cmin) round(cmax)])
ylabel({'Global','log_1_0 zmeso'})
title('Hist relationship')
lg  = legend({'CAN','CMCC','CNRM','GFDL','IPSL','UKESM','obsGLMM','obsSM','obsC1','obsC2'}); 
lg.Position(1:2) = [.75 .45];
lg.AutoUpdate = 'off';
text(-2,4.5,'A','FontWeight','Bold','FontSize',14)

subplot('Position',[0.45 0.765 0.25 0.2])
plot(Oy4,xs,'k','LineWidth',1.5); hold on;
plot(Sy4,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy4,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy4,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,ypred4,'--k','LineWidth',1.5); hold on;
for i=1:6
    plot(Hist.slop(i,4),log_diff(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0 2 -0.15 0.05])
ylabel('\Delta (log_1_0 zmeso)')
title('Future change')
text(0.02,0.0625,'B','FontWeight','Bold','FontSize',14)

% LC
subplot('Position',[0.1 0.525 0.25 0.2])
plot(log10(chl),hzoo1,'LineWidth',1.5); hold on;
ylim([-6 5])
xlim([round(cmin) round(cmax)])
ylabel({'LC','log_1_0 zmeso'})
text(-8.5,-35,'log_1_0 zmeso (mgC m^-^2)','Rotation',90)

subplot('Position',[0.45 0.525 0.25 0.2])
plot(x,x,':k'); hold on;
plot(Oy1,xs,'k','LineWidth',1.5); hold on;
plot(Sy1,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy1,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy1,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,ypred1,'--k','LineWidth',1.5); hold on;
for i=1:6
    plot(Hist.slop(i,1),log_diff(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0 2.5 -0.35 0.0])
ylabel('\Delta (log_1_0 zmeso)')

% HCSS
subplot('Position',[0.1 0.29 0.25 0.2])
plot(log10(chl),hzoo2,'LineWidth',1.5); hold on;
xlim([round(cmin) round(cmax)])
ylim([-6 5])
ylabel({'HCSS','log_1_0 zmeso'})

subplot('Position',[0.45 0.29 0.25 0.2])
plot(Oy2,xs,'k','LineWidth',1.5); hold on;
plot(Sy2,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy2,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy2,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,ypred2,'--k','LineWidth',1.5); hold on;
for i=1:6
    plot(Hist.slop(i,2),log_diff(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([-0.5 3.5 -0.2 0.0])
ylabel('\Delta (log_1_0 zmeso)')

% HCPS
subplot('Position',[0.1 0.05 0.25 0.2])
plot(log10(chl),hzoo3,'LineWidth',1.5); hold on;
ylim([-6 5])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')
ylabel({'HCPS','log_1_0 zmeso'})

subplot('Position',[0.45 0.05 0.25 0.2])
plot(Oy3,xs,'k','LineWidth',1.5); hold on;
plot(Sy3,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy3,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy3,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,ypred3,'--k','LineWidth',1.5); hold on;
for i=1:6
    plot(Hist.slop(i,3),log_diff(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([-0.5 3 -0.15 0.0])
xlabel('Historic relationship')
ylabel('\Delta (log_1_0 zmeso)')

print('-dpng',[figp 'Hist_lm_SSP585_delta_log10_EC_global_histbiomes_cope.png'])

%% No CAN + global
f4 = figure('Units','inches','Position',[1 1 9 12]);
% Global
%subplot(4,4,1)
subplot('Position',[0.1 0.765 0.25 0.2])
for i=2:10
    plot(log10(chl),hzoo4(i,:),'color',cb(i,:),'LineWidth',1.5); hold on;
end
ylim([-2 4])
xlim([round(cmin) round(cmax)])
ylabel({'Global','log_1_0 zmeso'})
title('Hist relationship')
lg  = legend({'CMCC','CNRM','GFDL','IPSL','UKESM','obsGLMM','obsSM','obsC1','obsC2'}); 
lg.Position(1:2) = [.75 .45];
lg.AutoUpdate = 'off';
text(-2,4.4,'A','FontWeight','Bold','FontSize',14)

subplot('Position',[0.45 0.765 0.25 0.2])
plot(Oy4,xs,'k','LineWidth',1.5); hold on;
plot(Sy4,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy4,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy4,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,noCypred4,'--k','LineWidth',1.5); hold on;
for i=2:6
    plot(Hist.slop(i,4),log_diff(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0.25 1.25 -0.14 0.0])
ylabel('\Delta (log_1_0 zmeso)')
title('Future change')
text(0.25,0.01,'B','FontWeight','Bold','FontSize',14)

% LC
subplot('Position',[0.1 0.525 0.25 0.2])
for i=2:10
    plot(log10(chl),hzoo1(i,:),'color',cb(i,:),'LineWidth',1.5); hold on;
end
ylim([-3 4])
xlim([round(cmin) round(cmax)])
ylabel({'LC','log_1_0 zmeso'})
text(-8.5,-35,'log_1_0 zmeso (mgC m^-^2)','Rotation',90)

subplot('Position',[0.45 0.525 0.25 0.2])
plot(x,x,':k'); hold on;
plot(Oy1,xs,'k','LineWidth',1.5); hold on;
plot(Sy1,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy1,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy1,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,noCypred1,'--k','LineWidth',1.5); hold on;
for i=2:6
    plot(Hist.slop(i,1),log_diff(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0 2.5 -0.35 0.0])
ylabel('\Delta (log_1_0 zmeso)')

% HCSS
subplot('Position',[0.1 0.29 0.25 0.2])
for i=2:10
    plot(log10(chl),hzoo2(i,:),'color',cb(i,:),'LineWidth',1.5); hold on;
end
xlim([round(cmin) round(cmax)])
ylim([-4 4])
ylabel({'HCSS','log_1_0 zmeso'})

subplot('Position',[0.45 0.29 0.25 0.2])
plot(Oy2,xs,'k','LineWidth',1.5); hold on;
plot(Sy2,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy2,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy2,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,noCypred2,'--k','LineWidth',1.5); hold on;
for i=2:6
    plot(Hist.slop(i,2),log_diff(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([0 2 -0.05 0.0])
ylabel('\Delta (log_1_0 zmeso)')

% HCPS
subplot('Position',[0.1 0.05 0.25 0.2])
for i=2:10
    plot(log10(chl),hzoo3(i,:),'color',cb(i,:),'LineWidth',1.5); hold on;
end
ylim([-1.5 4])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')
ylabel({'HCPS','log_1_0 zmeso'})

subplot('Position',[0.45 0.05 0.25 0.2])
plot(Oy3,xs,'k','LineWidth',1.5); hold on;
plot(Sy3,xs,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold on;
plot(CSy3,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(CMy3,xs,'color',[0.6 0.6 0.6],'LineWidth',1.5); hold on;
plot(xpred,noCypred3,'--k','LineWidth',1.5); hold on;
for i=2:6
    plot(Hist.slop(i,3),log_diff(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
axis([-0.5 1 -0.125 0.0])
xlabel('Historic relationship')
ylabel('\Delta (log_1_0 zmeso)')
print('-dpng',[figp 'Hist_lm_SSP585_delta_log10_EC_global_histbiomes_cope_noCAN.png'])


