% Calculate different skill metrics for each ESM
% log10 transform biomass
% Kelly Kearney skill calcs

clear all
close all

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

load([sfile 'skill_scores_hist_model_obsglm100_clim_log10_taylor.mat'])

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
%[rmsdA,itA]=sort(skill(4,2:7,1),'descend');  %CRMSD
rmsdA = skill(4,2:7,1);         %CRMSD
thetaA = acos(skill(2,2:7,1));  %R
rhoA = skill(6,2:7,1);          %norm std
% Winter
rmsdD = skill(4,2:7,2);
thetaD=acos(skill(2,2:7,2));
rhoD=skill(6,2:7,2);
% Spring
rmsdM = skill(4,2:7,3);
thetaM=acos(skill(2,2:7,3));
rhoM=skill(6,2:7,3);
% Summer
rmsdJ = skill(4,2:7,4);
thetaJ=acos(skill(2,2:7,4));
rhoJ=skill(6,2:7,4);
% Fall
rmsdS = skill(4,2:7,5);
thetaS=acos(skill(2,2:7,5));
rhoS=skill(6,2:7,5);

xc = 1; yc = 0;
th = linspace(0,2*pi,50);
r = 5:5:20;
rN = 1:20;
rP = 0.5:0.5:6;

%% 
%pp(ii)=polar(theta(ii),rho(ii));

% plotting points in non-polar?
% Peter Rochford:
%X = rho(2:end).*cos(theta(2:end));
%Y = rho(2:end).*sin(theta(2:end));

% built-in:
%rho   = STDs;
%theta = real(acos(CORs));
%pp(ii)=plot(rho(ii)*cos(theta(ii)),rho(ii)*sin(theta(ii)));

% Note that by definition the following relation must be true for all series i:
%RMSs(i) - sqrt(STDs(i).^2 + STDs(1)^2 - 2*STDs(i)*STDs(1).*CORs(i)) = 0

for j = 1:6
    for k = 1:5
        check(j,k) = skill(4,j+1,k) - sqrt(skill(6,j+1,k).^2 + skill(6,1,1)^2 ...
            - 2*skill(6,j+1,k)*skill(6,1,1).*skill(2,j+1,k));
    end
end

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
text(-0.15,2.5,'Standard deviation','HorizontalAlignment','center')
text(1.2,3.25,'Correlation coeff')
text(0.8,3.25,{'Root mean','square error'},'color',[0.5 0.5 0.5])
%title('log_1_0')
%text(2.2,4.7,'a','FontWeight','Bold','FontSize',14)

% t.Padding = 'compact';
% t.TileSpacing = 'compact';

print('-dpng',[figp 'Taylor_all_clim_together_obsglm_log10_ms_KK.png'])





