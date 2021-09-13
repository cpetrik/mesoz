% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zD','units');

cz = zD;
cz(cz(:)<0) = 0;
cunits = units;

clear zD units

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zD','units');

nz = zD;
nz(nz(:)<0) = 0;
nunits = units;

clear zD units

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zD','units');

uz = zD;
uz(uz(:)<0) = 0;
uunits = units;

clear zD units

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zD','units');

iz = zD;
iz(iz(:)<0) = 0;
iunits = units;

clear zD units

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso_200_climatol_1950_2014_gridded.mat'],...
    'zD','units');

gz = zD;
gz(gz(:)<0) = 0;
gunits = units;

clear zD units

%% COPEPOD
load('copepod-2012_cmass_DJF_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

% From m-3 to m-2 (integrate top 200m)
zoo_200 = zoo_g*200;

%% Vectorize, put in ABC order
model(:,1) = log10(cz(:));
model(:,2) = log10(nz(:));
model(:,3) = log10(gz(:));
model(:,4) = log10(iz(:));
model(:,5) = log10(uz(:));
obs = log10(zoo_200(:));

%% Skill metric = weighted sum of squares for multivariate

metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

skill=NaN*ones(14,5);
for j=1:5
    o=obs;
    p=model(:,j);
    
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

%% Results
simtex = {'CAN','CNRM','GFDL','IPSL','UK'};

% Bar graphs
figure(1)
subplot(3,2,1)
bar(skill(1,:),'k')
title('DJF Correlation coefficient')
set(gca,'XTickLabel',simtex)

subplot(3,2,3)
bar(skill(2,:),'k')
title('Root mean square error')
set(gca,'XTickLabel',simtex)

subplot(3,2,5)
bar(skill(6,:),'k')
title('Modeling efficiency')
set(gca,'XTickLabel',simtex)
ylim([-45 5])
print('-dpng','corr_rmse_mef_winter_clim_200.png')


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
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.35 0.35 0.35]}; %grey

tr=0;
rr=1;
figure(2)
subplot(1,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 5])
title('DJF')

subplot(1,2,2)
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
title('DJF')
legend([' ' simtex])
legend('location','east')
print('-dpng','Taylor_Joliff_winter_clim_200')

%% Taylor diagram using corr coeff calculated
theta2=acos(skill(1,:));

figure(3)
subplot(1,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 5])
title('DJF')

subplot(1,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 1.25])
title('DJF')
legend([' ' simtex])
legend('location','northeast')
print('-dpng','Taylor_winter_clim_200')


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

%%
nn = ~isnan(obs);
model = model(nn,:);
obs = obs(nn);

lat = lat_g(:);
lon = lon_g(:);
lat = lat(nn,:);
lon = lon(nn,:);

comb(:,1) = lat;
comb(:,2) = lon;
comb(:,3) = obs;
comb(:,4:8) = model;

obsmod = array2table(comb,'VariableNames',...
    {'Lat','Lon','obs','CAN','CNRM','GFDL','IPSL','UK'});

save('skill_winter_clim_200.mat','model','obs','metrics','skill','T','obsmod');
writetable(obsmod,'obs_mod_winter_clim_200.csv')

%% Scatter plots
%add corr line

% x=-10:0.1:5;
% 
% figure(4)
% subplot(2,3,1)
% scatter(obs,model(:,1)); hold on
% plot(x,x,'--k')
% axis([0 5 -10 4.5])
% 
% subplot(2,3,2)
% scatter(obs,model(:,2));hold on
% plot(x,x,'--k')
% axis([0 5 1.5 4.5])
% 
% subplot(2,3,3)
% scatter(obs,model(:,3));hold on
% plot(x,x,'--k')
% axis([0.5 5 1.5 4])
% 
% subplot(2,3,4)
% scatter(obs,model(:,4));hold on
% plot(x,x,'--k')
% axis([0 5 1 3.5])
% 
% subplot(2,3,5)
% scatter(obs,model(:,5));hold on
% plot(x,x,'--k')
% axis([0 5 0 4])








