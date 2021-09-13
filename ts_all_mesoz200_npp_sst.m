% CMIP6 output 
% 200m integrations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([cpath 'can_hist_ssp585_tsmeans_npp.mat']);

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([npath 'cnrm_hist_ssp585_tsmeans_npp.mat']);

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([upath 'ukesm_hist_ssp585_tsmeans_npp.mat']);

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([ipath 'ipsl_hist_ssp585_tsmeans_npp.mat']);

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([gpath 'gfdl_hist_ssp585_tsmeans_npp.mat']);

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/';
load([mpath 'cmcc_hist_ssp585_tsmeans_zmeso200_schl_sst.mat']);
load([mpath 'cmcc_hist_ssp585_tsmeans_npp.mat']);

%%
Czm = [HCzm(end-49:end) FCzm];
Mzm = [HMzm(end-49:end) FMzm];
Nzm = [HNzm(end-49:end) FNzm];
Gzm = [HGzm(end-49:end) FGzm];
Izm = [HIzm(end-49:end) FIzm];
Uzm = [HUzm(end-49:end) FUzm];

Cnpp = [HCnpp(end-49:end) FCnpp];
Mnpp = [HMnpp(end-49:end) FMnpp];
Nnpp = [HNnpp(end-49:end) FNnpp];
Gnpp = [HGnpp(end-49:end) FGnpp];
Inpp = [HInpp(end-49:end) FInpp];
Unpp = [HUnpp(end-49:end) FUnpp];

Csst = [HCsst(end-49:end) FCsst];
Msst = [HMsst(end-49:end) FMsst];
Nsst = [HNsst(end-49:end) FNsst];
Gsst = [HGsst(end-49:end) FGsst];
Isst = [HIsst(end-49:end) FIsst];
Usst = [HUsst(end-49:end) FUsst];

yr = 1965:2100;

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
Czm = Czm * 12.01 * 1e3;
Mzm = Mzm * 12.01 * 1e3;
Nzm = Nzm * 12.01 * 1e3;
Gzm = Gzm * 12.01 * 1e3;
Izm = Izm * 12.01 * 1e3;
Uzm = Uzm * 12.01 * 1e3;

%npp in molC/m2/s1 to g/m2/yr
%all models in molC: 12.01 g C in 1 mol C
Cnpp = Cnpp * 12.01 * 60*60*24*365;
Mnpp = Mnpp * 12.01 * 60*60*24*365;
Nnpp = Nnpp * 12.01 * 60*60*24*365;
Gnpp = Gnpp * 12.01 * 60*60*24*365;
Inpp = Inpp * 12.01 * 60*60*24*365;
Unpp = Unpp * 12.01 * 60*60*24*365;

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
cid = find(yr==1965);
mid = find(yr==1965);

%% 2 x 3 change order
figure('Units','inches','Position',[1 2 10 5]);

% Percent change from 1965
subplot(2,3,1)
plot(yr,100*(Czm-Czm(cid))/Czm(cid),'LineWidth',1.5); hold on;
plot(yr,100*(Mzm-Mzm(mid))/Mzm(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nzm-Nzm(yid))/Nzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gzm-Gzm(yid))/Gzm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Izm-Izm(yid))/Izm(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Uzm-Uzm(yid))/Uzm(yid),'LineWidth',1.5); hold on;
title('% \Delta Mesozoo')
ylabel('Percent change')
text(1965,12,'A')
xlim([1965 2100])

subplot(2,3,2)
plot(yr,100*(Cnpp-Cnpp(cid))/Cnpp(cid),'LineWidth',1.5); hold on;
plot(yr,100*(Mnpp-Mnpp(mid))/Mnpp(mid),'LineWidth',1.5); hold on;
plot(yr,100*(Nnpp-Nnpp(yid))/Nnpp(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Gnpp-Gnpp(yid))/Gnpp(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Inpp-Inpp(yid))/Inpp(yid),'LineWidth',1.5); hold on;
plot(yr,100*(Unpp-Unpp(yid))/Unpp(yid),'LineWidth',1.5); hold on;
title('% \Delta NPP')
ylabel('Percent change')
text(1965,32,'B')
xlim([1965 2100])

% Change from 1965
subplot(2,3,3)
plot(yr,Csst-Csst(cid),'LineWidth',1.5); hold on;
plot(yr,Msst-Msst(mid),'LineWidth',1.5); hold on;
plot(yr,Nsst-Nsst(yid),'LineWidth',1.5); hold on;
plot(yr,Gsst-Gsst(yid),'LineWidth',1.5); hold on;
plot(yr,Isst-Isst(yid),'LineWidth',1.5); hold on;
plot(yr,Usst-Usst(yid),'LineWidth',1.5); hold on;
title('\Delta SST')
ylabel('Change (^oC)')
text(1965,6.4,'C')
xlim([1965 2100])
legend({'CAN','CMCC','CNRM','GFDL','IPSL','UK'})
legend('location','northwest')

%Raw
subplot(2,3,4)
plot(yr,Czm,'LineWidth',1.5); hold on;
plot(yr,Mzm,'LineWidth',1.5); hold on;
plot(yr,Nzm,'LineWidth',1.5); hold on;
plot(yr,Gzm,'LineWidth',1.5); hold on;
plot(yr,Izm,'LineWidth',1.5); hold on;
plot(yr,Uzm,'LineWidth',1.5); hold on;
title('Mesozoo')
ylabel('Biomass (mgC m^-^2)')
text(1965,1080,'D')
xlim([1965 2100])

subplot(2,3,5)
plot(yr,Cnpp,'LineWidth',1.5); hold on;
plot(yr,Mnpp,'LineWidth',1.5); hold on;
plot(yr,Nnpp,'LineWidth',1.5); hold on;
plot(yr,Gnpp,'LineWidth',1.5); hold on;
plot(yr,Inpp,'LineWidth',1.5); hold on;
plot(yr,Unpp,'LineWidth',1.5); hold on;
title('NPP')
ylabel('Production (gC m^-^2 y^-^1)')
xlim([1965 2100])
text(1965,1.65e2,'E')

subplot(2,3,6)
plot(yr,Csst,'LineWidth',1.5); hold on;
plot(yr,Msst,'LineWidth',1.5); hold on;
plot(yr,Nsst,'LineWidth',1.5); hold on;
plot(yr,Gsst,'LineWidth',1.5); hold on;
plot(yr,Isst,'LineWidth',1.5); hold on;
plot(yr,Usst,'LineWidth',1.5); hold on;
title('SST')
ylabel('^oC')
text(1965,19,'F')
xlim([1965 2100])
print('-dpng',[ppath 'Pdiff_Raw_hist_ssp585_tsmeans_zmeso200_npp_sst.png'])




