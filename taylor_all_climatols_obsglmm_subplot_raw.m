% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

load('skill_model_obsglm_climatols.mat')
load('skill_scores_model_obsglm_clim_raw.mat')

%dim 1: metrics
%dim 2: model
%dim 3: season

figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

mtex = {'CAN','CNRM','GFDL','IPSL','UK'};
stext = {'All','Winter','Spring','Summer','Fall'};

%% Taylor diagram using Joliff corr coeff
cb=[34/255 136/255 51/255;...   %green
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

%% Annual
[rmsd,it]=sort(skill(10,:,1),'descend');  %tot RMSD
theta=acos(skill(7,:,1));                 %Taylor R
rho=skill(8,:,1);                         %Taylor norm std
mtex{6}='obs';

tr=0;
rr=1;
figure(1)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 12])
title('Annual climatology')
%%
subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 3])
title('Annual climatology')
legend([' ' mtex])
legend('location','eastoutside')
print('-dpng',[figp 'Taylor_Joliff_all_clim_200_obsglm_raw'])

%% Winter
[rmsd,it]=sort(skill(10,:,2),'descend');  %tot RMSD
theta=acos(skill(7,:,2));                 %Taylor R
rho=skill(8,:,2);                         %Taylor norm std

tr=0;
rr=1;
figure(2)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 18])
title('Winter')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 5])
title('Winter')
legend([' ' mtex])
legend('location','eastoutside')
print('-dpng',[figp 'Taylor_winter_clim_200_obsglm_raw'])

%% Spring
[rmsd,it]=sort(skill(10,:,3),'descend');  %tot RMSD
theta=acos(skill(7,:,3));                 %Taylor R
rho=skill(8,:,3); 

tr=0;
rr=1;
figure(3)
subplot(2,3,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 18])
title('Spring')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 6])
title('Spring')
legend([' ' mtex])
legend('location','eastoutside')
print('-dpng',[figp 'Taylor_spring_clim_200_obsglm_raw'])

%% Summer
[rmsd,it]=sort(skill(10,:,4),'descend');  %tot RMSD
theta=acos(skill(7,:,4));                 %Taylor R
rho=skill(8,:,4); 

tr=0;
rr=1;
figure(4)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 10])
title('Summer')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 3])
title('Summer')
legend([' ' mtex])
legend('location','eastoutside')
print('-dpng',[figp 'Taylor_summer_clim_200_obsglm_raw'])

%% Fall
[rmsd,it]=sort(skill(10,:,5),'descend');  %tot RMSD
theta=acos(skill(7,:,5));                 %Taylor R
rho=skill(8,:,5);  

tr=0;
rr=1;
figure(5)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 15])
title('Fall')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 4])
title('Fall')
legend([' ' mtex])
legend('location','eastoutside')
print('-dpng',[figp 'Taylor_fall_clim_200_obsglm_raw'])


%% All as subplots
% Annual
[rmsd,it]=sort(skill(10,:,1),'descend');  %tot RMSD
theta=acos(skill(7,:,1));                 %Taylor R
rho=skill(8,:,1);                         %Taylor norm std
mtex{6}='obs';
tr=0;
rr=1;

f10 = figure('Units','inches','Position',[1 3 7.5 10]);
subplot(5,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 12])
title('Annual climatology')

subplot(5,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 3])
title('Annual climatology')
%legend([' ' mtex])
%legend('location','eastoutside')


% Winter
[rmsd,it]=sort(skill(10,:,2),'descend');  %tot RMSD
theta=acos(skill(7,:,2));                 %Taylor R
rho=skill(8,:,2);                         %Taylor norm std

subplot(5,2,3)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 18])
title('Winter')

subplot(5,2,4)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 5])
title('Winter')
%legend([' ' mtex])
%legend('location','eastoutside')


% Spring
[rmsd,it]=sort(skill(10,:,3),'descend');  %tot RMSD
theta=acos(skill(7,:,3));                 %Taylor R
rho=skill(8,:,3); 

subplot(5,2,5)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 18])
title('Spring')
legend([' ' mtex])
legend('location','westoutside')

subplot(5,2,6)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 6])
title('Spring')
%legend([' ' mtex])
%legend('location','eastoutside')


% Summer
[rmsd,it]=sort(skill(10,:,4),'descend');  %tot RMSD
theta=acos(skill(7,:,4));                 %Taylor R
rho=skill(8,:,4); 

subplot(5,2,7)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 10])
title('Summer')

subplot(5,2,8)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 3])
title('Summer')
%legend([' ' mtex])
%legend('location','eastoutside')


% Fall
[rmsd,it]=sort(skill(10,:,5),'descend');  %tot RMSD
theta=acos(skill(7,:,5));                 %Taylor R
rho=skill(8,:,5);  

subplot(5,2,9)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 15])
title('Fall')

subplot(5,2,10)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 4])
title('Fall')
%legend([' ' mtex])
%legend('location','eastoutside')
print('-dpng',[figp 'Taylor_all_clim_200_obsglm_raw_subplots'])



%% Try combining on one
[rmsd,it]=sort(skill(10,:,5),'descend');  %tot RMSD
theta=acos(skill(7,:,5));                 %Taylor R
rho=skill(8,:,5);  

tr=0;
rr=1;
figure(5)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 15])
title('Fall')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cb(s,:),'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 4])
title('Fall')
legend([' ' mtex])
legend('location','eastoutside')
%print('-dpng',[figp 'Taylor_fall_clim_200_obsglm_raw'])






