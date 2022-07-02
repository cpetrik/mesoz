% Calculate different skill metrics for each ESM
% log10 transform biomass

clear all
close all

load('skill_scores_hist_model_obsglm100_clim_log10.mat')

%dim 1: metrics
%dim 2: model
%dim 3: season

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

mtex = {'obsGLMM','CAN','CMCC','CNRM','GFDL','IPSL','UK'};
stext = {'All','Winter','Spring','Summer','Fall'};
atex = {'obsGLMM','Annual','Winter','Spring','Summer','Fall',...
    'CAN','CMCC','CNRM','GFDL','IPSL','UK'};

rtl1 = (90 - [30,60,90])/100;
rtl2 = (100 - [30,60,90])/90;

%% Taylor diagram using Joliff corr coeff
% cb=[34/255 136/255 51/255;...   %green
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0.333 0.333 0.333];         %grey

cb=[34/255 136/255 51/255;...   %green
    153/255 153/255 51/255;...  %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];                     %black

set(groot,'defaultAxesColorOrder',cb);

%% RAW
tr=0;
rr=1;

thetaTickP = acos([0.9 0.6 0.3 0]);
thetaTickN = acos([0.9 0.6 0.3 0 -0.3 -0.6 -0.9]);

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

xc = 1; yc = 0;
th = linspace(0,2*pi,50);
r = 5:5:20;
rN = 1:20;
rP = 0.5:0.5:6;

%%
f1 = figure('Units','inches','Position',[1 1 8 6]);
%t = tiledlayout(2,1);
%nexttile

subplot('position',[0.1 0.15 0.8 0.8])
%plot obs
h2=polarplot(tr,rr,'k*'); hold on;
set(h2,'MarkerSize',10);
%plot all CAN to get shapes
for s=1
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(7,:),'MarkerSize',8);
    hD=polarplot(thetaD(s),rhoD(s),'d'); hold on;
    set(hD,'color','k','MarkerFaceColor',cb(7,:),'MarkerSize',8);
    hM=polarplot(thetaM(s),rhoM(s),'^'); hold on;
    set(hM,'color','k','MarkerFaceColor',cb(7,:),'MarkerSize',8);
    hJ=polarplot(thetaJ(s),rhoJ(s),'s'); hold on;
    set(hJ,'color','k','MarkerFaceColor',cb(7,:),'MarkerSize',8);
    hS=polarplot(thetaS(s),rhoS(s),'v'); hold on;
    set(hS,'color','k','MarkerFaceColor',cb(7,:),'MarkerSize',8);
end
%plot all annuals to get colors
for s=1:6
    hA=polarplot(thetaA(s),rhoA(s),'o'); hold on;
    set(hA,'color','k','MarkerFaceColor',cb(s,:),'MarkerSize',8);
end
%legend of colors and shapes
lgd = legend(atex,'Location','southoutside','NumColumns',2);
lgd.AutoUpdate = 'off';
%plot RMSE circles
for i = 1:7%length(rP)
    [x,y] = pol2cart(th,rP(i));
    [th1, r1] = cart2pol( x+xc, y+yc );
    hR=polarplot(th1, r1); hold on;
    %set(hR,'linestyle','--','color',[238/255 119/255 51/255]);
    set(hR,'linestyle','--','color',[0.7 0.7 0.7]);
    if (rem(i,2)==0)
        text(th1(22),r1(22),num2str(rP(i)),'color',[0.7 0.7 0.7]);
    end
end
%add all other data
for s=1:6
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
%relabel angles to match corr coeff
set(gca,'ThetaTick',(90/1.5708)*thetaTickN,'ThetaTickLabel',...
    {'0.9','0.6','0.3','0','-0.3','-0.6','-0.9'})
axis([0 180 0 2.5])
text(-0.15,2.0,'Standard deviation','HorizontalAlignment','center')
text(1.2,2.75,'Correlation coeff')
text(0.8,2.75,{'Root mean','square error'},'color',[0.5 0.5 0.5])
%title('log_1_0')
%text(2.2,4.7,'a','FontWeight','Bold','FontSize',14)

% t.Padding = 'compact';
% t.TileSpacing = 'compact';

print('-dpng',[figp 'Taylor_all_clim_together_obsglm_log10_ms.png'])




