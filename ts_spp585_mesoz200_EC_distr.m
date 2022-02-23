% CMIP6 output 
% 200m integrations
% plot time series of hist and ssp585
% plot mean of ensemble with and without EC reduction
% plot pdf of change in mesoz with and without EC reduction

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
load([mpath 'cmcc_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);

%%
Czm = [HCzm FCzm];
Mzm = [HMzm FMzm];
Nzm = [HNzm FNzm];
Gzm = [HGzm FGzm];
Izm = [HIzm FIzm];
Uzm = [HUzm FUzm];

Cyr = 1951:2100;
yr = 1950:2100;
Myr = 1965:2100;

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
Czm = Czm * 12.01 * 1e3;
Mzm = Mzm * 12.01 * 1e3;
Nzm = Nzm * 12.01 * 1e3;
Gzm = Gzm * 12.01 * 1e3;
Izm = Izm * 12.01 * 1e3;
Uzm = Uzm * 12.01 * 1e3;

%% Plots 
% cb=[34/255 136/255 51/255;...   %green
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0];                     %black
%     0.333 0.333 0.333];         %grey

cb=[34/255 136/255 51/255;...   %green
    %238/255 119/255 51/255;...  %orange
    %0/255 153/255 136/255;...   %teal
    153/255 153/255 51/255;...   %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];                     %black

set(groot,'defaultAxesColorOrder',cb);

%% Change from 1965
yid = find(yr==1965);
cid = find(Cyr==1965);
mid = find(Myr==1965);

pdiff = nan*ones(6,length(yr));

pdiff(1,2:end) = (Czm-Czm(cid))/Czm(cid);
pdiff(2,16:end) = (Mzm-Mzm(mid))/Mzm(mid);
pdiff(3,:) = (Nzm-Nzm(yid))/Nzm(yid);
pdiff(4,:) = (Gzm-Gzm(yid))/Gzm(yid);
pdiff(5,:) = (Izm-Izm(yid))/Izm(yid);
pdiff(6,:) = (Uzm-Uzm(yid))/Uzm(yid);

%% mean and stdev
Fpdiff = 100*pdiff(:,16:end);   %full ensemble
Epdiff = 100*pdiff(3:5,16:end); %best EC models

mF = mean(Fpdiff);
mE = mean(Epdiff);
sF = std(Fpdiff);
sE = std(Epdiff);

mo = 1965:2100;
X=[mo fliplr(mo)]; 
%create y values for out and then back 
ECI =[mE+sE, fliplr(mE-sE)]; 
FCI =[mF+sF, fliplr(mF-sF)]; 

% pdf contrained and not
Fpd = makedist('Normal','mu',mF(end),'sigma',sF(end));
Epd = makedist('Normal','mu',mE(end),'sigma',sE(end));
pc = -30:15;
Fpdf_norm = pdf(Fpd,pc);
Epdf_norm = pdf(Epd,pc);

%% save
save('ECdistr_fig_Pdiff_Raw_hist_ssp585_tsmeans_zmeso200.mat','yr','pdiff',...
    'Myr','mF','X','FCI','mE','ECI','pc','Fpdf_norm','Epdf_norm')

%% 2 x 3 change order
figure('Units','inches','Position',[1 1 16 4]);

% Percent change from 1965
subplot(1,3,1)
plot(yr,100*pdiff,'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
xlabel('Year')
text(1965,11,'A')
xlim([1965 2100])
legend({'CAN','CMCC','CNRM','GFDL','IPSL','UK'})
legend('location','southwest')

% Ensemble means of % change
subplot(1,3,2)
plot(Myr,mF,'-.k','LineWidth',1.5); hold on;
fill(X,FCI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
plot(Myr,mE,'k','LineWidth',1.5); hold on;
fill(X,ECI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
xlabel('Year')
text(1965,7,'B')
xlim([1965 2100])

% PDFs of % change
subplot(1,3,3)
plot(pc,Fpdf_norm,'-.k','LineWidth',2); hold on;
plot(pc,Epdf_norm,'k','LineWidth',2); hold on;
title('% \Delta Mesozoo')
xlabel('Percent change')
ylabel('Probability density')
text(-30,0.1,'C')
print('-dpng',[ppath 'Pdiff_Raw_hist_ssp585_tsmeans_zmeso200_ECdistr.png'])

%%
figure('Units','inches','Position',[1 1 10 4]);

% Percent change from 1965
subplot(1,2,1)
plot(Myr,mF,'-.k','LineWidth',1.5); hold on;
plot(Myr,mE,'k','LineWidth',1.5); hold on;
fill(X,FCI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
fill(X,ECI,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on;
ylabel('% \Delta Mesozoo')
xlabel('Year')
text(1965,6,'A','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
legend('full','constrained')
legend('location','southwest')

% PDFs of % change
subplot(1,2,2)
plot(pc,Fpdf_norm,'-.k','LineWidth',2); hold on;
plot(pc,Epdf_norm,'k','LineWidth',2); hold on;
xlabel('% \Delta Mesozoo')
ylabel('Probability density')
text(-30,0.073,'B','FontWeight','Bold','FontSize',14)
print('-dpng',[ppath 'Pdiff_hist_ssp585_zmeso200_ECdistr.png'])




