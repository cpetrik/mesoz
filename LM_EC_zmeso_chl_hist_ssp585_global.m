% Linear regressions of zmeso biomass with surf chl 
% Global, annual climatology
% Historic and SSP585

clear all
close all

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load('Hist_SSP585_coeffs_mod_chl_sst_biomes.mat','HistChlHCSS',...
    'HistChlHCPS','SspChlHCSS','SspChlHCPS');
load('Hist_SSP585_coeffs_mod_chl_sst_LC45.mat');
load('Hist_SSP585_coeffs_mod_chl_global_obsglm.mat')

Hist.slop = [HistChlLC45(:,2),HistChlHCSS(:,2),HistChlHCPS(:,2),HistChlGlobal(:,2)];
Hist.int = [HistChlLC45(:,1),HistChlHCSS(:,1),HistChlHCPS(:,1),HistChlGlobal(:,1)];

SSP.slop = [SspChlLC45(:,2),SspChlHCSS(:,2),SspChlHCPS(:,2),SspChlGlobal(:,2)];
SSP.int = [SspChlLC45(:,1),SspChlHCSS(:,1),SspChlHCPS(:,1),SspChlGlobal(:,1)];

%% predict mesozoo over chl values
% Models have a very large range
cmin = log(0.01);
cmax = log(25);
chl = logspace(cmin,cmax);

hzoo1 = Hist.int(:,1) + Hist.slop(:,1) .* log(chl);
szoo1 = SSP.int(:,1) + SSP.slop(:,1) .* log(chl);
hzoo2 = Hist.int(:,2) + Hist.slop(:,2) .* log(chl);
szoo2 = SSP.int(:,2) + SSP.slop(:,2) .* log(chl);
hzoo3 = Hist.int(:,3) + Hist.slop(:,3) .* log(chl);
szoo3 = SSP.int(:,3) + SSP.slop(:,3) .* log(chl);
hzoo4 = Hist.int(:,4) + Hist.slop(:,4) .* log(chl);
szoo4 = SSP.int(:,4) + SSP.slop(:,4) .* log(chl);

%% color same as Taylor diagrams
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

hmod = {'CAN','CNRM','GFDL','IPSL','UK','obs'};
smod = {'CAN','CNRM','GFDL','IPSL','UK'};

%% linear regression of hist vs. future 

mdl1 = fitlm(Hist.slop(1:5,1),SSP.slop(:,1));
ypred1 = predict(mdl1,Hist.slop(1:5,1));

mdl2 = fitlm(Hist.slop(1:5,2),SSP.slop(:,2));
ypred2 = predict(mdl2,Hist.slop(1:5,2));

mdl3 = fitlm(Hist.slop(1:5,3),SSP.slop(:,3));
ypred3 = predict(mdl3,Hist.slop(1:5,3));

mdl4 = fitlm(Hist.slop(1:5,4),SSP.slop(:,4));
ypred4 = predict(mdl4,Hist.slop(1:5,4));

noCmdl1 = fitlm(Hist.slop(2:5,1),SSP.slop(2:5,1));
noCypred1 = predict(noCmdl1,Hist.slop(2:5,1));

noCmdl2 = fitlm(Hist.slop(2:5,2),SSP.slop(2:5,2));
noCypred2 = predict(noCmdl2,Hist.slop(2:5,2));

noCmdl3 = fitlm(Hist.slop(2:5,3),SSP.slop(2:5,3));
noCypred3 = predict(noCmdl3,Hist.slop(2:5,3));

noCmdl4 = fitlm(Hist.slop(2:5,4),SSP.slop(2:5,4));
noCypred4 = predict(noCmdl4,Hist.slop(2:5,4));

test=mdl1.Coefficients;
tt=table2array(test);

%% Figures
x = (-5:0.5:4);
xs = 0:0.5:6;
y1 = Hist.slop(6,1)*ones(size(xs));
y2 = Hist.slop(6,2)*ones(size(xs));
y3 = Hist.slop(6,3)*ones(size(xs));
y4 = Hist.slop(6,4)*ones(size(xs));

%% 3 biomes + global
figure(1)
subplot(3,3,1)
plot(log(chl),hzoo4,'LineWidth',2)
ylim([-15 10])
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,2)
plot(log(chl),szoo4,'LineWidth',2)
ylim([-15 10])
xlim([round(cmin) round(cmax)])
title('SSP 585')
xlabel('log chl (mg m^-^3)')

subplot(3,3,3)
for i=1:5
    plot(Hist.slop(i,4),SSP.slop(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(y4,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(1:5,4),ypred4,'k'); hold on;
axis([0 2 0 2])
title('Future chl sensitivity')
xlabel('Historic chl sensitivity')

% No CAN + global
subplot(3,3,4)
for i=2:6
    plot(log(chl),hzoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-5 5])
xlim([round(cmin) round(cmax)])
title('Hist')
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')

subplot(3,3,5)
for i=2:5
    plot(log(chl),szoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-5 5])
xlim([round(cmin) round(cmax)])
title('SSP 585')
xlabel('log chl (mg m^-^3)')

subplot(3,3,6)
for i=2:5
    plot(Hist.slop(i,4),SSP.slop(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(y4,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(2:5,4),noCypred4,'k'); hold on;
axis([0.25 1.25 0.25 1.5])
title('Future chl sensitivity')
xlabel('Historic chl sensitivity')

% No CAN no UK + global
subplot(3,3,7)
for i=2:4
    plot(log(chl),hzoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
plot(log(chl),hzoo4(6,:),'color',cb(6,:),'LineWidth',2); hold on;
ylim([-3 2])
xlim([round(cmin) round(cmax)])
title('Hist')
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')

subplot(3,3,8)
for i=2:4
    plot(log(chl),szoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-3 2])
xlim([round(cmin) round(cmax)])
title('SSP 585')
xlabel('log chl (mg m^-^3)')

subplot(3,3,9)
for i=2:4
    plot(Hist.slop(i,4),SSP.slop(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(y4,xs,'color',cb(6,:)); hold on;
axis([0.35 0.75 0.35 0.75])
title('Future chl sensitivity')
xlabel('Historic chl sensitivity')
print('-dpng',[figp 'Hist_SSP585_lm_EC_zmeso_chl_global_only.png'])

