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

%%
Czm = [HCzm FCzm];
Nzm = [HNzm FNzm];
Gzm = [HGzm FGzm];
Izm = [HIzm FIzm];
Uzm = [HUzm FUzm];

Cchl = [HCchl FCchl];
Nchl = [HNchl FNchl];
Gchl = [HGchl FGchl];
Ichl = [HIchl FIchl];
Uchl = [HUchl FUchl];

Csst = [HCsst FCsst];
Nsst = [HNsst FNsst];
Gsst = [HGsst FGsst];
Isst = [HIsst FIsst];
Usst = [HUsst FUsst];

Cyr = 1951:2100;
yr = 1950:2100;

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
Czm = Czm * 12.01 * 1e3;
Nzm = Nzm * 12.01 * 1e3;
Gzm = Gzm * 12.01 * 1e3;
Izm = Izm * 12.01 * 1e3;
Uzm = Uzm * 12.01 * 1e3;

%chl in kg m-3', put in g m-3 (except CNRM & IPSL already that)
Cchl = Cchl * 1e3;
Gchl = Gchl * 1e3;
Uchl = Uchl * 1e3;

%% Plots 
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
    0 0 0];                     %black
%     0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

%%
% Change from 1965
yid = find(yr==1965);


%%
figure(1)
subplot(2,2,1)
plot(Cyr,Czm); hold on;
plot(yr,Nzm); hold on;
plot(yr,Gzm); hold on;
plot(yr,Izm); hold on;
plot(yr,Uzm); hold on;
title('mesozooplankton')
ylabel('biomass (mgC m^-^2)')
xlim([1951 2100])

subplot(2,2,3)
plot(Cyr,Cchl); hold on;
plot(yr,Nchl); hold on;
plot(yr,Gchl); hold on;
plot(yr,Ichl); hold on;
plot(yr,Uchl); hold on;
title('surf chl')
ylabel('concentration (g m-3)')
xlim([1951 2100])

subplot(2,2,4)
plot(Cyr,Csst); hold on;
plot(yr,Nsst); hold on;
plot(yr,Gsst); hold on;
plot(yr,Isst); hold on;
plot(yr,Usst); hold on;
title('sst')
ylabel('^oC')
xlim([1951 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')
print('-dpng',[ppath 'Raw_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%% log10
figure(2)
subplot(2,2,1)
plot(Cyr,log10(Czm)); hold on;
plot(yr,log10(Nzm)); hold on;
plot(yr,log10(Gzm)); hold on;
plot(yr,log10(Izm)); hold on;
plot(yr,log10(Uzm)); hold on;
title('mesozooplankton')
ylabel('biomass (log_1_0 mgC m^-^2)')
xlim([1951 2100])

subplot(2,2,3)
plot(Cyr,log10(Cchl)); hold on;
plot(yr,log10(Nchl)); hold on;
plot(yr,log10(Gchl)); hold on;
plot(yr,log10(Ichl)); hold on;
plot(yr,log10(Uchl)); hold on;
title('surf chl')
ylabel('concentration (log_1_0 g m-3)')
xlim([1951 2100])

subplot(2,2,4)
plot(Cyr,Csst); hold on;
plot(yr,Nsst); hold on;
plot(yr,Gsst); hold on;
plot(yr,Isst); hold on;
plot(yr,Usst); hold on;
title('sst')
ylabel('^oC')
xlim([1951 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')
print('-dpng',[ppath 'Log10_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%% Change from 1965
figure(3)
subplot(2,2,1)
plot(Cyr,Czm-Czm(yid)); hold on;
plot(yr,Nzm-Nzm(2)); hold on;
plot(yr,Gzm-Gzm(2)); hold on;
plot(yr,Izm-Izm(2)); hold on;
plot(yr,Uzm-Uzm(2)); hold on;
title('mesozooplankton')
ylabel('biomass (mgC m^-^2)')
xlim([1951 2100])

subplot(2,2,3)
plot(Cyr,Cchl-Cchl(yid)); hold on;
plot(yr,Nchl-Nchl(2)); hold on;
plot(yr,Gchl-Gchl(2)); hold on;
plot(yr,Ichl-Ichl(2)); hold on;
plot(yr,Uchl-Uchl(2)); hold on;
title('surf chl')
ylabel('concentration (g m-3)')
xlim([1951 2100])

subplot(2,2,4)
plot(Cyr,Csst-Csst(yid)); hold on;
plot(yr,Nsst-Nsst(2)); hold on;
plot(yr,Gsst-Gsst(2)); hold on;
plot(yr,Isst-Isst(2)); hold on;
plot(yr,Usst-Usst(2)); hold on;
title('sst')
ylabel('^oC')
xlim([1951 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')
print('-dpng',[ppath 'Rel_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])


%%
figure('Units','inches','Position',[1 4 10 5]);
subplot(2,3,1)
plot(Cyr,Czm); hold on;
plot(yr,Nzm); hold on;
plot(yr,Gzm); hold on;
plot(yr,Izm); hold on;
plot(yr,Uzm); hold on;
title('mesozoo')
ylabel('biomass (mgC m^-^2)')
xlim([1951 2100])

subplot(2,3,2)
plot(Cyr,Cchl*1e3); hold on;
plot(yr,Nchl*1e3); hold on;
plot(yr,Gchl*1e3); hold on;
plot(yr,Ichl*1e3); hold on;
plot(yr,Uchl*1e3); hold on;
title('surf chl')
ylabel('concentration (mg m-3)')
xlim([1951 2100])

subplot(2,3,3)
plot(Cyr,Csst); hold on;
plot(yr,Nsst); hold on;
plot(yr,Gsst); hold on;
plot(yr,Isst); hold on;
plot(yr,Usst); hold on;
title('sst')
ylabel('^oC')
xlim([1951 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

% Change from 1951
subplot(2,3,4)
plot(Cyr,Czm-Czm(yid)); hold on;
plot(yr,Nzm-Nzm(2)); hold on;
plot(yr,Gzm-Gzm(2)); hold on;
plot(yr,Izm-Izm(2)); hold on;
plot(yr,Uzm-Uzm(2)); hold on;
title('\Delta mesozoo')
ylabel('biomass (mgC m^-^2)')
xlim([1951 2100])

subplot(2,3,5)
plot(Cyr,(Cchl-Cchl(yid))*1e3); hold on;
plot(yr,(Nchl-Nchl(2))*1e3); hold on;
plot(yr,(Gchl-Gchl(2))*1e3); hold on;
plot(yr,(Ichl-Ichl(2))*1e3); hold on;
plot(yr,(Uchl-Uchl(2))*1e3); hold on;
title('\Delta surf chl')
ylabel('concentration (mg m-3)')
xlim([1951 2100])

subplot(2,3,6)
plot(Cyr,Csst-Csst(yid)); hold on;
plot(yr,Nsst-Nsst(2)); hold on;
plot(yr,Gsst-Gsst(2)); hold on;
plot(yr,Isst-Isst(2)); hold on;
plot(yr,Usst-Usst(2)); hold on;
title('\Delta sst')
ylabel('^oC')
xlim([1951 2100])
print('-dpng',[ppath 'Rel_Raw_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%% 3 x 3
figure('Units','inches','Position',[1 3 10 8]);
subplot(3,3,1)
plot(Cyr,Czm); hold on;
plot(yr,Nzm); hold on;
plot(yr,Gzm); hold on;
plot(yr,Izm); hold on;
plot(yr,Uzm); hold on;
title('mesozoo')
ylabel('mgC m^-^2')
xlim([1965 2100])

subplot(3,3,2)
plot(Cyr,Cchl*1e3); hold on;
plot(yr,Nchl*1e3); hold on;
plot(yr,Gchl*1e3); hold on;
plot(yr,Ichl*1e3); hold on;
plot(yr,Uchl*1e3); hold on;
title('surf chl')
ylabel('mg m-3')
xlim([1965 2100])

subplot(3,3,3)
plot(Cyr,Csst); hold on;
plot(yr,Nsst); hold on;
plot(yr,Gsst); hold on;
plot(yr,Isst); hold on;
plot(yr,Usst); hold on;
title('sst')
ylabel('^oC')
xlim([1965 2100])

% Change from 1965
subplot(3,3,4)
plot(Cyr,Czm-Czm(yid)); hold on;
plot(yr,Nzm-Nzm(yid)); hold on;
plot(yr,Gzm-Gzm(yid)); hold on;
plot(yr,Izm-Izm(yid)); hold on;
plot(yr,Uzm-Uzm(yid)); hold on;
title('\Delta mesozoo')
ylabel('mgC m^-^2')
xlim([1965 2100])

subplot(3,3,5)
plot(Cyr,(Cchl-Cchl(yid))*1e3); hold on;
plot(yr,(Nchl-Nchl(yid))*1e3); hold on;
plot(yr,(Gchl-Gchl(yid))*1e3); hold on;
plot(yr,(Ichl-Ichl(yid))*1e3); hold on;
plot(yr,(Uchl-Uchl(yid))*1e3); hold on;
title('\Delta surf chl')
ylabel('mg m-3')
xlim([1965 2100])

subplot(3,3,6)
plot(Cyr,Csst-Csst(yid)); hold on;
plot(yr,Nsst-Nsst(yid)); hold on;
plot(yr,Gsst-Gsst(yid)); hold on;
plot(yr,Isst-Isst(yid)); hold on;
plot(yr,Usst-Usst(yid)); hold on;
title('\Delta sst')
ylabel('^oC')
xlim([1965 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

% Percent change from 1965
subplot(3,3,7)
plot(Cyr,100*(Czm-Czm(yid))/Czm(yid)); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid)); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid)); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid)); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid)); hold on;
title('% \Delta mesozoo')
xlim([1965 2100])

%
subplot(3,3,8)
plot(Cyr,100*(Cchl-Cchl(yid))/Cchl(yid)); hold on;
plot(yr,100*(Nchl-Nchl(yid))/Nchl(yid)); hold on;
plot(yr,100*(Gchl-Gchl(yid))/Gchl(yid)); hold on;
plot(yr,100*(Ichl-Ichl(yid))/Ichl(yid)); hold on;
plot(yr,100*(Uchl-Uchl(yid))/Uchl(yid)); hold on;
title('% \Delta surf chl')
xlim([1965 2100])

% subplot(3,3,9)
% plot(Cyr,100*(Csst-Csst(yid))/Csst(yid)); hold on;
% plot(yr,100*(Nsst-Nsst(yid))/Nsst(yid)); hold on;
% plot(yr,100*(Gsst-Gsst(yid))/Gsst(yid)); hold on;
% plot(yr,100*(Isst-Isst(yid))/Isst(yid)); hold on;
% plot(yr,100*(Usst-Usst(yid))/Usst(yid)); hold on;
% title('% \Delta sst')
% xlim([1965 2100])
print('-dpng',[ppath 'Pdiff_Rel_Raw_hist_ssp585_tsmeans_zmeso200_schl_sst.png'])

%% 2 x 3
figure('Units','inches','Position',[1 2 10 5]);
subplot(2,3,4)
plot(Cyr,Czm,'LineWidth',1.5); hold on;
plot(yr,Nzm,'LineWidth',1.5); hold on;
plot(yr,Gzm,'LineWidth',1.5); hold on;
plot(yr,Izm,'LineWidth',1.5); hold on;
plot(yr,Uzm,'LineWidth',1.5); hold on;
title('Mesozoo')
ylabel('biomass (mgC m^-^2)')
xlabel('Year')
xlim([1965 2100])
text(1965,1080,'D')

subplot(2,3,5)
plot(Cyr,Cchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Nchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Gchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Ichl*1e3,'LineWidth',1.5); hold on;
plot(yr,Uchl*1e3,'LineWidth',1.5); hold on;
title('Surf Chl')
ylabel('concentration (mg m-3)')
xlabel('Year')
xlim([1965 2100])
text(1965,0.425,'E')

subplot(2,3,6)
plot(Cyr,Csst,'LineWidth',1.5); hold on;
plot(yr,Nsst,'LineWidth',1.5); hold on;
plot(yr,Gsst,'LineWidth',1.5); hold on;
plot(yr,Isst,'LineWidth',1.5); hold on;
plot(yr,Usst,'LineWidth',1.5); hold on;
title('SST')
ylabel('^oC')
xlabel('Year')
xlim([1965 2100])
text(1965,19,'F')

% Change from 1965
subplot(2,3,3)
plot(Cyr,Csst-Csst(yid),'LineWidth',1.5); hold on;
plot(yr,Nsst-Nsst(yid),'LineWidth',1.5); hold on;
plot(yr,Gsst-Gsst(yid),'LineWidth',1.5); hold on;
plot(yr,Isst-Isst(yid),'LineWidth',1.5); hold on;
plot(yr,Usst-Usst(yid),'LineWidth',1.5); hold on;
title('\Delta SST')
ylabel('change ^oC')
xlim([1965 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')
text(1965,6.4,'C')

% Percent change from 1965
subplot(2,3,1)
plot(Cyr,100*(Czm-Czm(yid))/Czm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid),'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('percent change')
xlim([1965 2100])
text(1965,6,'A')

subplot(2,3,2)
plot(Cyr,100*(Cchl-Cchl(yid))/Cchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Nchl-Nchl(yid))/Nchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gchl-Gchl(yid))/Gchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Ichl-Ichl(yid))/Ichl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uchl-Uchl(yid))/Uchl(yid),'LineWidth',1.5); hold on;
title('% \Delta Surf Chl')
ylabel('percent change')
xlim([1965 2100])
text(1965,45,'B')
print('-dpng',[ppath 'Pdiff_Raw_hist_ssp585_tsmeans_zmeso200_schl_sst_cb.png'])


%% 2 x 3 change order
figure('Units','inches','Position',[1 2 10 5]);
subplot(2,3,6)
plot(Cyr,Czm,'LineWidth',1.5); hold on;
plot(yr,Nzm,'LineWidth',1.5); hold on;
plot(yr,Gzm,'LineWidth',1.5); hold on;
plot(yr,Izm,'LineWidth',1.5); hold on;
plot(yr,Uzm,'LineWidth',1.5); hold on;
title('Mesozoo')
ylabel('Biomass (mgC m^-^2)')
xlim([1965 2100])

subplot(2,3,5)
plot(Cyr,Cchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Nchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Gchl*1e3,'LineWidth',1.5); hold on;
plot(yr,Ichl*1e3,'LineWidth',1.5); hold on;
plot(yr,Uchl*1e3,'LineWidth',1.5); hold on;
title('Surf Chl')
ylabel('Concentration (mg m-3)')
xlim([1965 2100])
text(1960,1100,'D')

subplot(2,3,4)
plot(Cyr,Csst,'LineWidth',1.5); hold on;
plot(yr,Nsst,'LineWidth',1.5); hold on;
plot(yr,Gsst,'LineWidth',1.5); hold on;
plot(yr,Isst,'LineWidth',1.5); hold on;
plot(yr,Usst,'LineWidth',1.5); hold on;
title('SST')
ylabel('^oC')
xlim([1965 2100])

% Change from 1965
subplot(2,3,1)
plot(Cyr,Csst-Csst(yid),'LineWidth',1.5); hold on;
plot(yr,Nsst-Nsst(yid),'LineWidth',1.5); hold on;
plot(yr,Gsst-Gsst(yid),'LineWidth',1.5); hold on;
plot(yr,Isst-Isst(yid),'LineWidth',1.5); hold on;
plot(yr,Usst-Usst(yid),'LineWidth',1.5); hold on;
title('\Delta SST')
ylabel('Change (^oC)')
xlim([1965 2100])
legend({'CAN','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

% Percent change from 1965
subplot(2,3,3)
plot(Cyr,100*(Czm-Czm(yid))/Czm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid),'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
xlim([1965 2100])

%
subplot(2,3,2)
plot(Cyr,100*(Cchl-Cchl(yid))/Cchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Nchl-Nchl(yid))/Nchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gchl-Gchl(yid))/Gchl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Ichl-Ichl(yid))/Ichl(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uchl-Uchl(yid))/Uchl(yid),'LineWidth',1.5); hold on;
title('% \Delta Surf Chl')
ylabel('Percent change')
xlim([1965 2100])
print('-dpng',[ppath 'Pdiff_Raw_hist_ssp585_tsmeans_zmeso200_schl_sst_cb_v2.png'])




