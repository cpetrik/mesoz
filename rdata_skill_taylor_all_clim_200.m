% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON','units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

czmo_all = zmo_all;
czmo_DJF = zmo_DJF;
czmo_JJA = zmo_JJA;
czmo_MAM = zmo_MAM;
czmo_SON = zmo_SON;

cunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON','units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

nzmo_all = zmo_all;
nzmo_DJF = zmo_DJF;
nzmo_JJA = zmo_JJA;
nzmo_MAM = zmo_MAM;
nzmo_SON = zmo_SON;

nunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON')%,'units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

uzmo_all = zmo_all;
uzmo_DJF = zmo_DJF;
uzmo_JJA = zmo_JJA;
uzmo_MAM = zmo_MAM;
uzmo_SON = zmo_SON;

%uunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'zmo_all','zmo_DJF','zmo_JJA','zmo_MAM','zmo_SON');%,'units');

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

izmo_all = zmo_all;
izmo_DJF = zmo_DJF;
izmo_JJA = zmo_JJA;
izmo_MAM = zmo_MAM;
izmo_SON = zmo_SON;

%iunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat']);

zmo_all(zmo_all(:)<0) = 0;
zmo_DJF(zmo_DJF(:)<0) = 0;
zmo_JJA(zmo_JJA(:)<0) = 0;
zmo_MAM(zmo_MAM(:)<0) = 0;
zmo_SON(zmo_SON(:)<0) = 0;

gzmo_all = zmo_all;
gzmo_DJF = zmo_DJF;
gzmo_JJA = zmo_JJA;
gzmo_MAM = zmo_MAM;
gzmo_SON = zmo_SON;

%gunits = units;

clear zmo_all zmo_DJF zmo_JJA zmo_MAM zmo_SON units

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_mesoz.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj] = size(gzmo_all);
glmz_mo = reshape(glmobsmesoz,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);

%% climatologies
ozmo_DJF = nanmean(glmz_mo(:,:,[1 2 12]),3);
ozmo_MAM = nanmean(glmz_mo(:,:,3:5),3);
ozmo_JJA = nanmean(glmz_mo(:,:,6:8),3);
ozmo_SON = nanmean(glmz_mo(:,:,9:11),3);
ozmo_all = nanmean(glmz_mo,3);

%% check orientations
figure(1)
pcolor(czmo_all); shading flat;
title('CAN')

figure(2)
pcolor(nzmo_all); shading flat;
title('CNRM')

figure(3)
pcolor(izmo_all); shading flat;
title('IPSL')

figure(4)
pcolor(gzmo_all); shading flat;
title('GFDL')

figure(5)
pcolor(uzmo_all); shading flat;
title('UK')

figure(6)
pcolor(ozmo_all); shading flat;
title('Obs')

%% Flip
%UK, CNRM, CAN
czmo_all = fliplr(czmo_all);
czmo_DJF = fliplr(czmo_DJF);
czmo_JJA = fliplr(czmo_JJA);
czmo_MAM = fliplr(czmo_MAM);
czmo_SON = fliplr(czmo_SON);

nzmo_all = fliplr(nzmo_all);
nzmo_DJF = fliplr(nzmo_DJF);
nzmo_JJA = fliplr(nzmo_JJA);
nzmo_MAM = fliplr(nzmo_MAM);
nzmo_SON = fliplr(nzmo_SON);

uzmo_all = fliplr(uzmo_all);
uzmo_DJF = fliplr(uzmo_DJF);
uzmo_JJA = fliplr(uzmo_JJA);
uzmo_MAM = fliplr(uzmo_MAM);
uzmo_SON = fliplr(uzmo_SON);

%% Vectorize, put in ABC order
% convert to same units
%models mol C/m2 -> g C/m2
%obsglm mg C/m2 -> g C/m2
mod_all(:,1) = (czmo_all(:)) * 12.01;
mod_all(:,2) = (nzmo_all(:)) * 12.01;
mod_all(:,3) = (gzmo_all(:)) * 12.01;
mod_all(:,4) = (izmo_all(:)) * 12.01;
mod_all(:,5) = (uzmo_all(:)) * 12.01;
mod_all(:,6) = (ozmo_all(:)) * 1e-3;

mod_DJF(:,1) = (czmo_DJF(:)) * 12.01;
mod_DJF(:,2) = (nzmo_DJF(:)) * 12.01;
mod_DJF(:,3) = (gzmo_DJF(:)) * 12.01;
mod_DJF(:,4) = (izmo_DJF(:)) * 12.01;
mod_DJF(:,5) = (uzmo_DJF(:)) * 12.01;
mod_DJF(:,6) = (ozmo_DJF(:)) * 1e-3;

mod_JJA(:,1) = (czmo_JJA(:)) * 12.01;
mod_JJA(:,2) = (nzmo_JJA(:)) * 12.01;
mod_JJA(:,3) = (gzmo_JJA(:)) * 12.01;
mod_JJA(:,4) = (izmo_JJA(:)) * 12.01;
mod_JJA(:,5) = (uzmo_JJA(:)) * 12.01;
mod_JJA(:,6) = (ozmo_JJA(:)) * 1e-3;

mod_MAM(:,1) = (czmo_MAM(:)) * 12.01;
mod_MAM(:,2) = (nzmo_MAM(:)) * 12.01;
mod_MAM(:,3) = (gzmo_MAM(:)) * 12.01;
mod_MAM(:,4) = (izmo_MAM(:)) * 12.01;
mod_MAM(:,5) = (uzmo_MAM(:)) * 12.01;
mod_MAM(:,6) = (ozmo_MAM(:)) * 1e-3;

mod_SON(:,1) = (czmo_SON(:)) * 12.01;
mod_SON(:,2) = (nzmo_SON(:)) * 12.01;
mod_SON(:,3) = (gzmo_SON(:)) * 12.01;
mod_SON(:,4) = (izmo_SON(:)) * 12.01;
mod_SON(:,5) = (uzmo_SON(:)) * 12.01;
mod_SON(:,6) = (ozmo_SON(:)) * 1e-3;

%% all clim
nn = ~isnan(mod_all(:,6));
mod_all = mod_all(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3:8) = mod_all;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','CAN','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(obsmod,'skill_model_obsglm_all_clim_200.csv')

%% Winter
nn = ~isnan(mod_DJF(:,6));
mod_DJF = mod_DJF(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

dcomb(:,1) = lat;
dcomb(:,2) = lon;
dcomb(:,3:8) = mod_DJF;

omDJF = array2table(dcomb,'VariableNames',...
    {'Lat','Lon','CAN','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(omDJF,'skill_model_obsglm_DJF_clim_200.csv')

%% JJA clim
nn = ~isnan(mod_JJA(:,6));
mod_JJA = mod_JJA(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

jcomb(:,1) = lat;
jcomb(:,2) = lon;
jcomb(:,3:8) = mod_JJA;

omJJA = array2table(jcomb,'VariableNames',...
    {'Lat','Lon','CAN','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(omJJA,'skill_model_obsglm_JJA_clim_200.csv')

%% MAM clim
nn = ~isnan(mod_MAM(:,6));
mod_MAM = mod_MAM(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

mcomb(:,1) = lat;
mcomb(:,2) = lon;
mcomb(:,3:8) = mod_MAM;

omMAM = array2table(mcomb,'VariableNames',...
    {'Lat','Lon','CAN','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(omMAM,'skill_model_obsglm_MAM_clim_200.csv')

%% SON clim
nn = ~isnan(mod_SON(:,6));
mod_SON = mod_SON(nn,:);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

scomb(:,1) = lat;
scomb(:,2) = lon;
scomb(:,3:8) = mod_SON;

omSON = array2table(scomb,'VariableNames',...
    {'Lat','Lon','CAN','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(omSON,'skill_model_obsglm_SON_clim_200.csv')

%% Skill metric = weighted sum of squares for multivariate

metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

skill=NaN*ones(14,5);
for j=1:5
    % compare logs
    o=log(mod_all(:,6) + eps);
    p=log(mod_all(:,j) + eps);
    
    % Only compare model & obs with non-NaNs
    p(isnan(o))=[];
    o(isnan(o))=[];
    o(isnan(p))=[];
    p(isnan(p))=[];
    
    n = length(o);
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    skill(1,j) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    skill(2,j) = sqrt(num/n);
    
    % reliability index  %%What about predicted zero values?
    %q1=nansum((log(o./p)).^2);
    q1=nansum((o-p).^2);
    skill(3,j) = exp(sqrt(q1/n));
    
    % average error
    skill(4,j) = nansum(p-o) / n;
    
    % average absolute error
    skill(5,j) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    skill(6,j) = (num1-num2)/num1;
    
    % Taylor R
    num=nansum((o-omean).*(p-pmean));
    skill(7,j) = num/(n*osig*psig);
    
    % Taylor normalized std
    skill(8,j) = psig/osig;
    
    % unbiased root mean square difference
    % sign tells difference between model and obs std
    q1=nansum(((p-pmean)-(o-omean)).^2);
    if (psig>=osig)
        skill(9,j) = sqrt(q1/n);
    else
        skill(9,j) = -1*sqrt(q1/n);
    end
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    skill(10,j) = sqrt(q1/n);
    
    % normalized bias
    skill(11,j) = (pmean(1,1)-omean(1,1));%./osig;
    
    % Joliff et al 2009 S1
    R=skill(7,j);
    s=skill(8,j);
    num=2*(1+R);
    den=(s + (1/s)).^2;
    skill(12,j) = 1 - (num/den);
    
    % Joliff et al 2009 S2
    num=(1+R).^4;
    den=4*(s + (1/s)).^2;
    skill(13,j) = 1 - (num/den);
    
    % Joliff et al 2009 S3
    q1=exp(-((s-1).^2) / 0.18);
    q2=(1+R)/2;
    skill(14,j) = 1 - (q1*q2);
    
end

%% save 
simtex = {'CAN','CNRM','GFDL','IPSL','UK'};

save('skill_model_obsglm_clim_zmeso200.mat','comb','dcomb','jcomb','mcomb',...
    'scomb','skill','metrics','simtex')

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

% Bar graphs
figure(1)
subplot(3,2,1)
bar(skill(1,:),'k')
title('Annual Correlation coefficient')
set(gca,'XTickLabel',simtex)

subplot(3,2,3)
bar(skill(2,:),'k')
title('Root mean square error')
set(gca,'XTickLabel',simtex)

subplot(3,2,5)
bar(skill(6,:),'k')
title('Modeling efficiency')
set(gca,'XTickLabel',simtex)
ylim([-15 5])
print('-dpng',[figp 'corr_rmse_mef_all_clim_200_obsglm_raw.png'])


%% Taylor diagram using Joliff corr coeff
% cm={[1 0.5 0],...   %orange
%     [0.5 0.5 0],... %tan/army
%     [0 0.7 0],...   %g
%     [0 1 1],...     %c
%     [0 0 0.75],...  %b
%     [0.5 0 1],...   %purple
%     [1 0 1],...     %m
%     [1 0 0],...     %r
%     [0.5 0 0],...   %maroon
%     [0.35 0.35 0.35]}; %grey


[rmsd,it]=sort(skill(10,:),'descend');  %tot RMSD
theta=acos(skill(7,:));                 %Taylor R
rho=skill(8,:);                         %Taylor norm std
simtex{6}='obs';

%%
cm={[0 0.7 0],...   %g
    [0 0 0.65],...  %b
    [0.7 0 1],...   %purple
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.35 0.35 0.35]}; %grey

tr=0;
rr=1;
figure(2)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 12])
title('Annual climatology')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 3])
title('Annual climatology')
legend([' ' simtex])
legend('location','east')

subplot(2,2,3)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 1.25])
title('Annual climatology')
legend([' ' simtex])
legend('location','east')
print('-dpng',[figp 'Taylor_Joliff_all_clim_200_obsglm_raw'])

%% Taylor diagram using corr coeff calculated
theta2=acos(skill(1,:));

figure(3)
subplot(2,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 12])
title('Annual climatology')

subplot(2,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 180 0 3])
title('Annual climatology')

subplot(2,2,3)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 1.5])
title('Annual climatology')
legend([' ' simtex])
legend('location','south')
print('-dpng',[figp 'Taylor_all_clim_200_obsglm_raw'])


%% Best
one=[1;3;6;7;8];
zer=[2;4;5;9;10;11];
high=[12;13;14];
for z=1:14
    Z1=sum(z==one);
    Z2=sum(z==zer);
    Z3=sum(z==high);
    if (Z1>0)
        best{z,1}=simtex{find(min(abs(1-(skill(z,:))))==abs(1-(skill(z,:))))};
    elseif (Z2>0)
        best{z,1}=simtex{find(min(abs(0-(skill(z,:))))==abs(0-(skill(z,:))))};
    elseif (Z3>0)
        best{z,1}=simtex{find(max(skill(z,:))==skill(z,:))};
    end
end

T=table(metrics,best);

%% Scatter plots
%add corr line

x=-3:0.1:2;
mod_all = log(mod_all + eps);
obs=mod_all(:,6);

%%
figure(4)
subplot(2,3,1)
scatter(obs,mod_all(:,1)); hold on
plot(x,x,'--k')
%axis([0 5 -10 4.5])

subplot(2,3,2)
scatter(obs,mod_all(:,2));hold on
plot(x,x,'--k')
%axis([0 5 1.5 4.5])

subplot(2,3,3)
scatter(obs,mod_all(:,3));hold on
plot(x,x,'--k')
%axis([0.5 5 1.5 4])

subplot(2,3,4)
scatter(obs,mod_all(:,4));hold on
plot(x,x,'--k')
%axis([0 5 1 3.5])

subplot(2,3,5)
scatter(obs,mod_all(:,5));hold on
plot(x,x,'--k')
%axis([0 5 0 4])








