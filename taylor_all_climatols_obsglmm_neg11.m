% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

load('skill_model_obsglm_climatols.mat')
load('skill_scores_model_obsglm_clim_neg11.mat')

%dim 1: metrics
%dim 2: model
%dim 3: season

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

mtex = {'obs','CAN','CNRM','GFDL','IPSL','UK'};
stext = {'All','Winter','Spring','Summer','Fall'};
atex = {'obs','CAN','CNRM','GFDL','IPSL','UK','All','Winter','Spring',...
    'Summer','Fall'};

%% Taylor diagram using Joliff corr coeff
cb=[34/255 136/255 51/255;...   %green
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

%% -1 to 1
tr=0;
rr=1;

% Annual
[rmsdA,itA]=sort(skill(10,:,1),'descend');  %tot RMSD
thetaA=acos(skill(7,:,1));                 %Taylor R
rhoA=skill(8,:,1);                         %Taylor norm std
% Winter
[rmsdD,itD]=sort(skill(10,:,2),'descend');  
thetaD=acos(skill(7,:,2));                 
rhoD=skill(8,:,2);   
% Spring
[rmsdM,itM]=sort(skill(10,:,3),'descend');  
thetaM=acos(skill(7,:,3));                 
rhoM=skill(8,:,3); 
% Summer
[rmsdJ,itJ]=sort(skill(10,:,4),'descend');  
thetaJ=acos(skill(7,:,4));                 
rhoJ=skill(8,:,4); 
% Fall
[rmsdS,itS]=sort(skill(10,:,5),'descend');  
thetaS=acos(skill(7,:,5));                 
rhoS=skill(8,:,5);  

%%
f1 = figure('Units','inches','Position',[1 1 10 5]);
subplot(1,2,1)
h2=polarplot(tr,rr,'k*'); hold on;
set(h2,'MarkerSize',10);
for s=1:5
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
end
for s=1
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hD=polarplot(thetaD(s),rhoD(s),'d'); hold on;
    set(hD,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hM=polarplot(thetaM(s),rhoM(s),'^'); hold on;
    set(hM,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hJ=polarplot(thetaJ(s),rhoJ(s),'s'); hold on;
    set(hJ,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hS=polarplot(thetaS(s),rhoS(s),'v'); hold on;
    set(hS,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
end
lgd = legend(atex,'Location','southoutside','NumColumns',2);
lgd.AutoUpdate = 'off';
for s=2:5
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hD=polarplot(thetaD(s),rhoD(s),'d'); hold on;
    set(hD,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hM=polarplot(thetaM(s),rhoM(s),'^'); hold on;
    set(hM,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hJ=polarplot(thetaJ(s),rhoJ(s),'s'); hold on;
    set(hJ,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hS=polarplot(thetaS(s),rhoS(s),'v'); hold on;
    set(hS,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
end
axis([0 180 0 5])

subplot(1,2,2)
h2=polarplot(tr,rr,'k*'); hold on;
set(h2,'MarkerSize',10);
for s=1:5
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hD=polarplot(thetaD(s),rhoD(s),'d'); hold on;
    set(hD,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hM=polarplot(thetaM(s),rhoM(s),'^'); hold on;
    set(hM,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hJ=polarplot(thetaJ(s),rhoJ(s),'s'); hold on;
    set(hJ,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
    hS=polarplot(thetaS(s),rhoS(s),'v'); hold on;
    set(hS,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
end
%legend
axis([0 90 0 4])

print('-dpng',[figp 'Taylor_all_clim_together_obsglm_neg11'])





