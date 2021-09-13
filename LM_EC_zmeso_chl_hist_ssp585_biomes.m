% Linear regressions of zmeso biomass with surf chl 
% Tropical regions, annual climatology
% Historic and SSP585

clear all
close all

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load('Hist_SSP585_coeffs_mod_chl_sst_biomes.mat');

Hist.slop = [HistChlLC(:,2), HistChlHCSS(:,2), HistChlHCPS(:,2)];
Hist.int = [HistChlLC(:,1), HistChlHCSS(:,1), HistChlHCPS(:,1)];

SSP.slop = [SspChlLC(:,2), SspChlHCSS(:,2), SspChlHCPS(:,2)];
SSP.int = [SspChlLC(:,1), SspChlHCSS(:,1), SspChlHCPS(:,1)];

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

%% figures
figure(1)
subplot(3,2,1)
plot(log(chl),hzoo1,'LineWidth',2)
legend(hmod)
legend('location','northwest')
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('LC Hist')

subplot(3,2,2)
plot(log(chl),szoo1,'LineWidth',2)
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,2,3)
plot(log(chl),hzoo2,'LineWidth',2)
% legend(SSP.mod)
% legend('location','northwest')
xlim([round(cmin) round(cmax)])
% ylim([0 5])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('HCSS Hist')

subplot(3,2,4)
plot(log(chl),szoo2,'LineWidth',2)
xlim([round(cmin) round(cmax)])
% ylim([0 5])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,2,5)
plot(log(chl),hzoo3,'LineWidth',2)
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('HCPS Hist')

subplot(3,2,6)
plot(log(chl),szoo3,'LineWidth',2)
xlim([round(cmin) round(cmax)])
% ylim([0 5])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')
print('-dpng',[figp 'Hist_SSP585_lm_zmeso_chl_biomes.png'])

%% linear regression of hist vs. future 

mdl1 = fitlm(Hist.slop(1:5,1),SSP.slop(:,1));
ypred1 = predict(mdl1,Hist.slop(1:5,1));

mdl2 = fitlm(Hist.slop(1:5,2),SSP.slop(:,2));
ypred2 = predict(mdl2,Hist.slop(1:5,2));

mdl3 = fitlm(Hist.slop(1:5,3),SSP.slop(:,3));
ypred3 = predict(mdl3,Hist.slop(1:5,3));

noCmdl1 = fitlm(Hist.slop(2:5,1),SSP.slop(2:5,1));
noCypred1 = predict(noCmdl1,Hist.slop(2:5,1));

noCmdl2 = fitlm(Hist.slop(2:5,2),SSP.slop(2:5,2));
noCypred2 = predict(noCmdl2,Hist.slop(2:5,2));

noCmdl3 = fitlm(Hist.slop(2:5,3),SSP.slop(2:5,3));
noCypred3 = predict(noCmdl3,Hist.slop(2:5,3));


test=mdl1.Coefficients;
tt=table2array(test);

%%
x = (-5:0.5:4);
xs = 0:0.5:6;
y1 = Hist.slop(6,1)*ones(size(xs));
y2 = Hist.slop(6,2)*ones(size(xs));
y3 = Hist.slop(6,3)*ones(size(xs));


%% 3 biomes
figure(1)
subplot(2,2,1)
for i=1:5
    plot(Hist.slop(i,1),SSP.slop(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(Hist.slop(1:5,1),ypred1,'k'); hold on;
legend(smod)
legend('location','southeast')
axis([0 3 0 2])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('LC Zmeso vs. chl')

subplot(2,2,2)
plot(x,x,':k'); hold on;
for i=1:5
    plot(Hist.slop(i,2),SSP.slop(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(Hist.slop(1:5,2),ypred2,'k'); hold on;
% legend(smod)
% legend('location','northwest')
axis([0 4 0 6])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCSS Zmeso vs. chl')

subplot(2,2,3)
plot(x,x,':k'); hold on;
for i=1:5
    plot(Hist.slop(i,3),SSP.slop(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(Hist.slop(1:5,3),ypred3,'k'); hold on;
% legend(smod)
% legend('location','northwest')
axis([-1 3 0 5])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCPS Zmeso vs. chl')
print('-dpng',[figp 'Hist_SSP585_EC_zmeso_chl_biomes.png'])

%% 3 biomes
figure(2)
subplot(3,3,1)
plot(log(chl),hzoo1,'LineWidth',2)
% legend(hmod)
% legend('location','northwest')
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,2)
plot(log(chl),szoo1,'LineWidth',2)
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,3)
for i=1:5
    plot(Hist.slop(i,1),SSP.slop(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(y1,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(1:5,1),ypred1,'k'); hold on;
% legend(smod)
% legend('location','southeast')
axis([0 3 0 2])
% xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('LC')

subplot(3,3,4)
plot(log(chl),hzoo2,'LineWidth',2)
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,5)
plot(log(chl),szoo2,'LineWidth',2)
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,6)
plot(x,x,':k'); hold on;
for i=1:5
    plot(Hist.slop(i,2),SSP.slop(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y2,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(1:5,2),ypred2,'k'); hold on;
axis([0 4 0 6])
% xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCSS')

subplot(3,3,7)
plot(log(chl),hzoo3,'LineWidth',2)
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,8)
plot(log(chl),szoo3,'LineWidth',2)
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,9)
plot(x,x,':k'); hold on;
for i=1:5
    plot(Hist.slop(i,3),SSP.slop(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y3,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(1:5,3),ypred3,'k'); hold on;
axis([-1 3 0 5])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCPS')

print('-dpng',[figp 'Hist_SSP585_lm_EC_zmeso_chl_biomes.png'])

%% No CAN
figure(4)
subplot(3,3,1)
for i=2:6
    plot(log(chl),hzoo1(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
% legend(hmod)
% legend('location','northwest')
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,2)
for i=2:5
    plot(log(chl),szoo1(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,3)
for i=2:5
    plot(Hist.slop(i,1),SSP.slop(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(x,x,':k'); hold on;
plot(y1,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(2:5,1),noCypred1,'k'); hold on;
% legend(smod)
% legend('location','southeast')
axis([0 3 0 2])
% xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('LC')

subplot(3,3,4)
for i=2:6
    plot(log(chl),hzoo2(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,5)
for i=2:5
    plot(log(chl),szoo2(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
xlim([round(cmin) round(cmax)])
% xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,6)
plot(x,x,':k'); hold on;
for i=2:5
    plot(Hist.slop(i,2),SSP.slop(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y2,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(2:5,2),noCypred2,'k'); hold on;
axis([0 2 0 2])
% xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCSS')

subplot(3,3,7)
for i=2:6
    plot(log(chl),hzoo3(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
ylabel('log zmeso (mgC m^-^2)')
title('Hist')

subplot(3,3,8)
for i=2:5
    plot(log(chl),szoo3(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
xlim([round(cmin) round(cmax)])
xlabel('log chl (mg m^-^3)')
% ylabel('log zmeso (mgC m^-^2)')
title('SSP 585')

subplot(3,3,9)
plot(x,x,':k'); hold on;
for i=2:5
    plot(Hist.slop(i,3),SSP.slop(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y3,xs,'color',cb(6,:)); hold on;
plot(Hist.slop(2:5,3),noCypred3,'k'); hold on;
axis([-1 1 0 2])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('HCPS')

print('-dpng',[figp 'Hist_SSP585_lm_EC_zmeso_chl_biomes_noCAN.png'])

