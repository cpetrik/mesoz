% CMIP6 output 
% 200m integrations

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

Cchl = [HCchl FCchl];
Mchl = [HMchl FMchl];
Nchl = [HNchl FNchl];
Gchl = [HGchl FGchl];
Ichl = [HIchl FIchl];
Uchl = [HUchl FUchl];

Csst = [HCsst FCsst];
Msst = [HMsst FMsst];
Nsst = [HNsst FNsst];
Gsst = [HGsst FGsst];
Isst = [HIsst FIsst];
Usst = [HUsst FUsst];

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

%chl in kg m-3', put in g m-3 (except CNRM & IPSL already that)
Cchl = Cchl * 1e3;
Mchl = Mchl * 1e3;
Gchl = Gchl * 1e3;
Uchl = Uchl * 1e3;

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

%%
% Change from 1965
yid = find(yr==1965);
cid = find(Cyr==1965);
mid = find(Myr==1965);

%% 2 x 3 change order
figure('Units','inches','Position',[1 2 10 5]);

% Percent change from 1965
subplot(2,3,1)
plot(Cyr,100*(Czm-Czm(cid))/Czm(cid),'LineWidth',1.5); hold on;
plot(Myr,100*(Mzm-Mzm(mid))/Mzm(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid),'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
text(1965,12,'a','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

subplot(2,3,2)
plot(Cyr,100*(Cchl-Cchl(cid))/Cchl(cid),'LineWidth',1.5); hold on;
plot(Myr,100*(Mchl-Mchl(mid))/Mchl(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nchl-Nchl(yid))/Nchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gchl-Gchl(yid))/Gchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Ichl-Ichl(yid))/Ichl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uchl-Uchl(yid))/Uchl(yid),'LineWidth',1.5); hold on;
title('% \Delta Surf Chl')
ylabel('Percent change')
text(1965,45,'b','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

% Change from 1965
subplot(2,3,3)
plot(Cyr,Csst-Csst(cid),'LineWidth',1.5); hold on;
plot(Myr,Msst-Msst(mid),'LineWidth',1.5); hold on;
plot(yr,Nsst-Nsst(yid),'LineWidth',1.5); hold on;
plot(yr,Gsst-Gsst(yid),'LineWidth',1.5); hold on;
plot(yr,Isst-Isst(yid),'LineWidth',1.5); hold on;
plot(yr,Usst-Usst(yid),'LineWidth',1.5); hold on;
title('\Delta SST')
ylabel('Change (^oC)')
text(1965,6.4,'c','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
legend({'CAN','CMCC','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

%Raw
subplot(2,3,4)
plot(Cyr,Czm,'LineWidth',1.5); hold on;
plot(Myr,Mzm,'LineWidth',1.5); hold on;
plot(yr,Nzm,'LineWidth',1.5); hold on;
plot(yr,Gzm,'LineWidth',1.5); hold on;
plot(yr,Izm,'LineWidth',1.5); hold on;
plot(yr,Uzm,'LineWidth',1.5); hold on;
title('Mesozoo')
ylabel('Biomass (mgC m^-^2)')
text(1965,1080,'d','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

subplot(2,3,5)
plot(Cyr,Cchl*1e3,'LineWidth',1.5); hold on;
plot(Myr,Mchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Nchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Gchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Ichl*1e3,'LineWidth',1.5); hold on;
plot(yr,Uchl*1e3,'LineWidth',1.5); hold on;
title('Surf Chl')
ylabel('Concentration (mg m-3)')
xlim([1965 2100])
text(1965,0.425,'e','FontWeight','Bold','FontSize',14)

subplot(2,3,6)
plot(Cyr,Csst,'LineWidth',1.5); hold on;
plot(Myr,Msst,'LineWidth',1.5); hold on;
plot(yr,Nsst,'LineWidth',1.5); hold on;
plot(yr,Gsst,'LineWidth',1.5); hold on;
plot(yr,Isst,'LineWidth',1.5); hold on;
plot(yr,Usst,'LineWidth',1.5); hold on;
title('SST')
ylabel('^oC')
text(1965,19,'f','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
print('-dpng',[ppath 'Pdiff_Raw_hist_ssp585_tsmeans_zmeso200_schl_sst_cb_ms.png'])

%% 2 x 3 log10 concentration
figure('Units','inches','Position',[1 2 10 5]);

% Percent change from 1965
subplot(2,3,1)
plot(Cyr,100*(Czm-Czm(cid))/Czm(cid),'LineWidth',1.5); hold on;
plot(Myr,100*(Mzm-Mzm(mid))/Mzm(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid),'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
text(1965,12,'a','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

subplot(2,3,2)
plot(Cyr,100*(Cchl-Cchl(cid))/Cchl(cid),'LineWidth',1.5); hold on;
plot(Myr,100*(Mchl-Mchl(mid))/Mchl(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nchl-Nchl(yid))/Nchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gchl-Gchl(yid))/Gchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Ichl-Ichl(yid))/Ichl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uchl-Uchl(yid))/Uchl(yid),'LineWidth',1.5); hold on;
title('% \Delta Surf Chl')
ylabel('Percent change')
text(1965,45,'b','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

% Change from 1965
subplot(2,3,3)
plot(Cyr,Csst-Csst(cid),'LineWidth',1.5); hold on;
plot(Myr,Msst-Msst(mid),'LineWidth',1.5); hold on;
plot(yr,Nsst-Nsst(yid),'LineWidth',1.5); hold on;
plot(yr,Gsst-Gsst(yid),'LineWidth',1.5); hold on;
plot(yr,Isst-Isst(yid),'LineWidth',1.5); hold on;
plot(yr,Usst-Usst(yid),'LineWidth',1.5); hold on;
title('\Delta SST')
ylabel('Change (^oC)')
text(1965,6.4,'c','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
legend({'CAN','CMCC','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

%Raw
subplot(2,3,4)
plot(Cyr,log10(Czm),'LineWidth',1.5); hold on;
plot(Myr,log10(Mzm),'LineWidth',1.5); hold on;
plot(yr,log10(Nzm),'LineWidth',1.5); hold on;
plot(yr,log10(Gzm),'LineWidth',1.5); hold on;
plot(yr,log10(Izm),'LineWidth',1.5); hold on;
plot(yr,log10(Uzm),'LineWidth',1.5); hold on;
title('Mesozoo')
ylabel('Biomass (log_1_0 mgC m^-^2)')
text(1965,1080,'d','FontWeight','Bold','FontSize',14)
xlim([1965 2100])

subplot(2,3,5)
plot(Cyr,log10(Cchl*1e3),'LineWidth',1.5); hold on;
plot(Myr,log10(Mchl*1e3),'LineWidth',1.5); hold on;
plot(yr,log10(Nchl*1e3),'LineWidth',1.5); hold on;
plot(yr,log10(Gchl*1e3),'LineWidth',1.5); hold on;
plot(yr,log10(Ichl*1e3),'LineWidth',1.5); hold on;
plot(yr,log10(Uchl*1e3),'LineWidth',1.5); hold on;
title('Surf Chl')
ylabel('Concentration (log_1_0 mg m-3)')
xlim([1965 2100])
text(1965,0.425,'e','FontWeight','Bold','FontSize',14)

subplot(2,3,6)
plot(Cyr,Csst,'LineWidth',1.5); hold on;
plot(Myr,Msst,'LineWidth',1.5); hold on;
plot(yr,Nsst,'LineWidth',1.5); hold on;
plot(yr,Gsst,'LineWidth',1.5); hold on;
plot(yr,Isst,'LineWidth',1.5); hold on;
plot(yr,Usst,'LineWidth',1.5); hold on;
title('SST')
ylabel('^oC')
text(1965,19,'f','FontWeight','Bold','FontSize',14)
xlim([1965 2100])
print('-dpng',[ppath 'Pdiff_log10_hist_ssp585_tsmeans_zmeso200_schl_sst_cb_ms.png'])




