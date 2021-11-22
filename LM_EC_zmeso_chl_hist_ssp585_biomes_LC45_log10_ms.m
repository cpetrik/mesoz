% Linear regressions of zmeso biomass with surf chl 
% Annual climatology, biomes and global
% Historic and SSP585
% GLMM100

clear all
close all

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('Hist_SSP585_coeffs_log10_chl_global_biomes.mat');
% load('Hist_CIs_biomes_obsglm100.mat');
% load('Hist_std_err_biomes_obsglm100.mat');
load('Hist_CIs_std_err_biomes_stromberg.mat');

%%
Hist.slop = [HistB1,HistB2,HistB3,HistBG];
Hist.int = [HistA1,HistA2,HistA3,HistAG];

SSP.slop = [SSP585B1(1:6),SSP585B2(1:6),SSP585B3(1:6),SSP585BG(1:6)];
SSP.int = [SSP585A1(1:6),SSP585A2(1:6),SSP585A3(1:6),SSP585AG(1:6)];

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

%% replace hzoo4 stromberg
hzoo1(7,:) = sCI1(:,1)';
hzoo2(7,:) = sCI2(:,1)';
hzoo3(7,:) = sCI3(:,1)';
hzoo4(7,:) = sCIG(:,1)';

% replace hzoo4 w/glmm100
% hzoo1(7,:) = CI1(:,1)';
% hzoo2(7,:) = CI2(:,1)';
% hzoo3(7,:) = CI3(:,1)';
% hzoo4(7,:) = CIG(:,1)';

%% Confidence intervals as shadded region
%create continuous x value array for plotting
mo = log10(chl);
X=[mo fliplr(mo)]; 
%create y values for out and then back
hzoo1CI =[sCI1(:,3); flipud(sCI1(:,2))]; 
hzoo2CI =[sCI2(:,3); flipud(sCI2(:,2))]; 
hzoo3CI =[sCI3(:,3); flipud(sCI3(:,2))]; 
hzoo4CI =[sCIG(:,3); flipud(sCIG(:,2))]; 

%create y values for out and then back glmm100
% hzoo1CI =[CI1(:,3); flipud(CI1(:,2))]; 
% hzoo2CI =[CI2(:,3); flipud(CI2(:,2))]; 
% hzoo3CI =[CI3(:,3); flipud(CI3(:,2))]; 
% hzoo4CI =[CIG(:,3); flipud(CIG(:,2))]; 

%% color same as Taylor diagrams
cb=[34/255 136/255 51/255;...   %green
    153/255 153/255 51/255;...  %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];         %black
%     0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

%hmod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obs'};
hmod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','SM'};
smod = {'CAN','CMCC','CNRM','GFDL','IPSL','UK'};

%% linear regression of hist vs. future 

mdl1 = fitlm(Hist.slop(1:6,1),SSP.slop(:,1));
ypred1 = predict(mdl1,Hist.slop(1:6,1));

mdl2 = fitlm(Hist.slop(1:6,2),SSP.slop(:,2));
ypred2 = predict(mdl2,Hist.slop(1:6,2));

mdl3 = fitlm(Hist.slop(1:6,3),SSP.slop(:,3));
ypred3 = predict(mdl3,Hist.slop(1:6,3));

mdl4 = fitlm(Hist.slop(1:6,4),SSP.slop(:,4));
ypred4 = predict(mdl4,Hist.slop(1:6,4));

noCmdl1 = fitlm(Hist.slop(2:6,1),SSP.slop(2:6,1));
noCypred1 = predict(noCmdl1,Hist.slop(2:6,1));

noCmdl2 = fitlm(Hist.slop(2:6,2),SSP.slop(2:6,2));
noCypred2 = predict(noCmdl2,Hist.slop(2:6,2));

noCmdl3 = fitlm(Hist.slop(2:6,3),SSP.slop(2:6,3));
noCypred3 = predict(noCmdl3,Hist.slop(2:6,3));

noCmdl4 = fitlm(Hist.slop(2:6,4),SSP.slop(2:6,4));
noCypred4 = predict(noCmdl4,Hist.slop(2:6,4));

test=mdl1.Coefficients;
tt=table2array(test);

%% Figures
x = (-1:0.5:6);
xs = -1:0.5:6;

%glmm100
% y1 = bOG(2)*ones(size(xs));
% y2 = bOG(3)*ones(size(xs));
% y3 = bOG(4)*ones(size(xs));
% y4 = bOG(1)*ones(size(xs));
% 
% l1 = bOG(2)*ones(size(xs)) - 2*bOGse(2);
% l2 = bOG(3)*ones(size(xs)) - 2*bOGse(3);
% l3 = bOG(4)*ones(size(xs)) - 2*bOGse(4);
% l4 = bOG(1)*ones(size(xs)) - 2*bOGse(1);
% 
% u1 = bOG(2)*ones(size(xs)) + 2*bOGse(2);
% u2 = bOG(3)*ones(size(xs)) + 2*bOGse(3);
% u3 = bOG(4)*ones(size(xs)) + 2*bOGse(4);
% u4 = bOG(1)*ones(size(xs)) + 2*bOGse(1);

%stromberg
y1 = bSM(2)*ones(size(xs));
y2 = bSM(3)*ones(size(xs));
y3 = bSM(4)*ones(size(xs));
y4 = bSM(1)*ones(size(xs));

l1 = bSM(2)*ones(size(xs)) - 2*bSMse(2);
l2 = bSM(3)*ones(size(xs)) - 2*bSMse(3);
l3 = bSM(4)*ones(size(xs)) - 2*bSMse(4);
l4 = bSM(1)*ones(size(xs)) - 2*bSMse(1);

u1 = bSM(2)*ones(size(xs)) + 2*bSMse(2);
u2 = bSM(3)*ones(size(xs)) + 2*bSMse(3);
u3 = bSM(4)*ones(size(xs)) + 2*bSMse(4);
u4 = bSM(1)*ones(size(xs)) + 2*bSMse(1);

%% 3 biomes + global
f1 = figure('Units','inches','Position',[1 1 9 8]);
subplot(4,4,1)
plot(log10(chl),hzoo4,'LineWidth',2); hold on;
ylim([-6 3])
xlim([round(cmin) round(cmax)])
ylabel({'Global','log_1_0 zmeso'})
title('Hist')
lg  = legend({'CAN','CMCC','CNRM','GFDL','IPSL','UKESM','SM'}); 
lg.Position(1:2) = [.725 .45];
lg.AutoUpdate = 'off';
fill(X,hzoo4CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25) %plot filled area

subplot(4,4,2)
plot(log10(chl),szoo4,'LineWidth',2)
ylim([-6 3])
xlim([round(cmin) round(cmax)])
title('SSP 585')

subplot(4,4,3)
plot(x,x,'--k'); hold on;
for i=1:6
    plot(Hist.slop(i,4),SSP.slop(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y4,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l4,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u4,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(1:6,4),ypred4,'k'); hold on;
axis([0.2 2 0.2 2])
title('Future chl sensitivity')

subplot(4,4,5)
plot(log10(chl),hzoo1,'LineWidth',2); hold on;
fill(X,hzoo1CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
ylim([-6 5])
xlim([round(cmin) round(cmax)])
% xlabel('log10 chl (mg m^-^3)')
ylabel({'LC','log_1_0 zmeso'})
text(-8.5,-35,'log_1_0 zmeso (mgC m^-^2)','Rotation',90)

subplot(4,4,6)
plot(log10(chl),szoo1,'LineWidth',2)
ylim([-6 5])
xlim([round(cmin) round(cmax)])
% xlabel('log_1_0 chl (mg m^-^3)')
% ylabel('log_1_0 zmeso (mgC m^-^2)')
%title('SSP 585')

subplot(4,4,7)
plot(x,x,'--k'); hold on;
for i=1:6
    plot(Hist.slop(i,1),SSP.slop(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y1,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l1,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u1,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(1:6,1),ypred1,'k'); hold on;
% legend(smod)
% legend('location','southeast')
axis([0 2.5 0 2.5])

subplot(4,4,9)
plot(log10(chl),hzoo2,'LineWidth',2); hold on;
fill(X,hzoo2CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
xlim([round(cmin) round(cmax)])
ylim([-10 7])
% xlabel('log_1_0 chl (mg m^-^3)')
%ylabel('log_1_0 zmeso (mgC m^-^2)')
ylabel({'HCSS','log_1_0 zmeso'})

subplot(4,4,10)
plot(log10(chl),szoo2,'LineWidth',2)
xlim([round(cmin) round(cmax)])
ylim([-10 7])

subplot(4,4,11)
plot(x,x,'--k'); hold on;
for i=1:6
    plot(Hist.slop(i,2),SSP.slop(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y2,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l2,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u2,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(1:6,2),ypred2,'k'); hold on;
axis([0 4 0 6])
% xlabel('Historic chl sensitivity')
%ylabel('Future chl sensitivity')
%title('HCSS')

subplot(4,4,13)
plot(log10(chl),hzoo3,'LineWidth',2); hold on;
fill(X,hzoo1CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
ylim([-10 5])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')
ylabel({'HCPS','log_1_0 zmeso'})

subplot(4,4,14)
plot(log10(chl),szoo3,'LineWidth',2)
ylim([-10 5])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')
% ylabel('log_1_0 zmeso (mgC m^-^2)')
%title('SSP 585')

subplot(4,4,15)
plot(x,x,'--k'); hold on;
for i=1:6
    plot(Hist.slop(i,3),SSP.slop(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y3,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l3,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u3,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(1:6,3),ypred3,'k'); hold on;
axis([-0.5 3 0 5])
xlabel('Historic chl sensitivity')
%ylabel('Future chl sensitivity')
%title('HCPS')

print('-dpng',[figp 'Hist_SSP585_lm_EC_zmeso_chl_global_biomes_LC45_stromberg_CIs.png'])

%% No CAN + global
f4 = figure('Units','inches','Position',[1 1 9 8]);
subplot(4,4,1)
for i=2:7
    plot(log10(chl),hzoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-2 2])
lg  = legend({'CMCC','CNRM','GFDL','IPSL','UKESM','SM'}); 
lg.Position(1:2) = [.71 .45];
lg.AutoUpdate = 'off';
fill(X,hzoo4CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
xlim([round(cmin) round(cmax)])
title('Hist')
ylabel({'Global','log_1_0 zmeso'})

subplot(4,4,2)
for i=2:6
    plot(log10(chl),szoo4(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-2 2])
xlim([round(cmin) round(cmax)])
title('SSP 585')

subplot(4,4,3)
plot(x,x,'--k'); hold on;
for i=2:6
    plot(Hist.slop(i,4),SSP.slop(i,4),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y4,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l4,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u4,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(2:6,4),noCypred4,'k'); hold on;
axis([0.25 1.25 0.25 1.5])
title('Future chl sensitivity')

%LC
subplot(4,4,5)
for i=2:7
    plot(log10(chl),hzoo1(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
fill(X,hzoo1CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
ylim([-3 4])
xlim([round(cmin) round(cmax)])
ylabel({'LC','log_1_0 zmeso'})
text(-8.5,-30,'log_1_0 zmeso (mgC m^-^2)','Rotation',90)

subplot(4,4,6)
for i=2:6
    plot(log10(chl),szoo1(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-3 4])
xlim([round(cmin) round(cmax)])

subplot(4,4,7)
plot(x,x,'--k'); hold on;
for i=2:6
    plot(Hist.slop(i,1),SSP.slop(i,1),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y1,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l1,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u1,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(2:6,1),noCypred1,'k'); hold on;
axis([0 2.5 0 2.5])

%HCSS
subplot(4,4,9)
for i=2:7
    plot(log10(chl),hzoo2(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
fill(X,hzoo2CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
ylim([-4 3])
xlim([round(cmin) round(cmax)])
ylabel({'HCSS','log_1_0 zmeso'})

subplot(4,4,10)
for i=2:6
    plot(log10(chl),szoo2(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-4 3])
xlim([round(cmin) round(cmax)])

subplot(4,4,11)
plot(x,x,'--k'); hold on;
for i=2:6
    plot(Hist.slop(i,2),SSP.slop(i,2),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y2,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l2,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u2,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(2:6,2),noCypred2,'k'); hold on;
axis([-0.2 2 -0.2 1.4])

%HCPS
subplot(4,4,13)
for i=2:7
    plot(log10(chl),hzoo3(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
fill(X,hzoo3CI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25)
ylim([-1.5 1])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')
ylabel({'HCPS','log_1_0 zmeso'})

subplot(4,4,14)
for i=2:6
    plot(log10(chl),szoo3(i,:),'color',cb(i,:),'LineWidth',2); hold on;
end
ylim([-1.5 1])
xlim([round(cmin) round(cmax)])
xlabel('log_1_0 chl (mg m^-^3)')

subplot(4,4,15)
plot(x,x,'--k'); hold on;
for i=2:6
    plot(Hist.slop(i,3),SSP.slop(i,3),'.','color',cb(i,:),'MarkerSize',20); hold on;
end
plot(y3,xs,'color',[0.35 0.35 0.35]); hold on;
plot(l3,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(u3,xs,':','color',[0.35 0.35 0.35]); hold on;
plot(Hist.slop(2:6,3),noCypred3,'k'); hold on;
axis([-0.5 0.75 0 0.75])
xlabel('Historic chl sensitivity')
print('-dpng',[figp 'Hist_SSP585_lm_EC_zmeso_chl_global_biomes_LC45_stromberg_CIs_noCAN.png'])


