% CMIP6 output 
% 200m integrations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('Hist_ts_chl_clim_biomes.mat')
load('Hist_ts_sst_clim_biomes.mat')
load('Hist_ts_mesoz_clim_biomes.mat')

%% Plots 
% cm=[0 0.7 0;...   %g
%     0 0 0.75;...  %b
%     0.5 0 1;...   %purple
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
%     0.35 0.35 0.35]; %grey
cb=[51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0 0 0];                     %black
%     0.333 0.333 0.333];         %grey

set(groot,'defaultAxesColorOrder',cb);

%%
mo = 1:12;

close all

figure(1)
% MESOZ ------------------------------
subplot(4,3,1)
plot(mo,Nz,'LineWidth',1.5); hold on;
plot(mo,Gz,'LineWidth',1.5); hold on;
plot(mo,Iz,'LineWidth',1.5); hold on;
plot(mo,Uz,'LineWidth',1.5); hold on;
plot(mo,Oz,'LineWidth',1.5); hold on;
title('mesozoo (mgC m^-^2)')
ylabel('Global')
xlim([1 12])

subplot(4,3,4)
plot(mo,Nzlc,'LineWidth',1.5); hold on;
plot(mo,Gzlc,'LineWidth',1.5); hold on;
plot(mo,Izlc,'LineWidth',1.5); hold on;
plot(mo,Uzlc,'LineWidth',1.5); hold on;
plot(mo,Ozlc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,7)
plot(mo,Nzss,'LineWidth',1.5); hold on;
plot(mo,Gzss,'LineWidth',1.5); hold on;
plot(mo,Izss,'LineWidth',1.5); hold on;
plot(mo,Uzss,'LineWidth',1.5); hold on;
plot(mo,Ozss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,10)
plot(mo,Nzps,'LineWidth',1.5); hold on;
plot(mo,Gzps,'LineWidth',1.5); hold on;
plot(mo,Izps,'LineWidth',1.5); hold on;
plot(mo,Uzps,'LineWidth',1.5); hold on;
plot(mo,Ozps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])

% CHL --------------------------------
subplot(4,3,2)
plot(mo,Nc,'LineWidth',1.5); hold on;
plot(mo,Gc,'LineWidth',1.5); hold on;
plot(mo,Ic,'LineWidth',1.5); hold on;
plot(mo,Uc,'LineWidth',1.5); hold on;
plot(mo,Oc,'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
ylabel('Global')
xlim([1 12])

subplot(4,3,5)
plot(mo,Nclc,'LineWidth',1.5); hold on;
plot(mo,Gclc,'LineWidth',1.5); hold on;
plot(mo,Iclc,'LineWidth',1.5); hold on;
plot(mo,Uclc,'LineWidth',1.5); hold on;
plot(mo,Oclc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,8)
plot(mo,Ncss,'LineWidth',1.5); hold on;
plot(mo,Gcss,'LineWidth',1.5); hold on;
plot(mo,Icss,'LineWidth',1.5); hold on;
plot(mo,Ucss,'LineWidth',1.5); hold on;
plot(mo,Ocss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,11)
plot(mo,Ncps,'LineWidth',1.5); hold on;
plot(mo,Gcps,'LineWidth',1.5); hold on;
plot(mo,Icps,'LineWidth',1.5); hold on;
plot(mo,Ucps,'LineWidth',1.5); hold on;
plot(mo,Ocps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlabel('NH month')
xlim([1 12])

% SST --------------------------
subplot(4,3,3)
plot(mo,Nt,'LineWidth',1.5); hold on;
plot(mo,Gt,'LineWidth',1.5); hold on;
plot(mo,It,'LineWidth',1.5); hold on;
plot(mo,Ut,'LineWidth',1.5); hold on;
plot(mo,Ot,'LineWidth',1.5); hold on;
title('SST (^oC)')
ylabel('Global')
xlim([1 12])

subplot(4,3,6)
plot(mo,Ntlc,'LineWidth',1.5); hold on;
plot(mo,Gtlc,'LineWidth',1.5); hold on;
plot(mo,Itlc,'LineWidth',1.5); hold on;
plot(mo,Utlc,'LineWidth',1.5); hold on;
plot(mo,Otlc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

subplot(4,3,9)
plot(mo,Ntss,'LineWidth',1.5); hold on;
plot(mo,Gtss,'LineWidth',1.5); hold on;
plot(mo,Itss,'LineWidth',1.5); hold on;
plot(mo,Utss,'LineWidth',1.5); hold on;
plot(mo,Otss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

subplot(4,3,12)
plot(mo,Ntps,'LineWidth',1.5); hold on;
plot(mo,Gtps,'LineWidth',1.5); hold on;
plot(mo,Itps,'LineWidth',1.5); hold on;
plot(mo,Utps,'LineWidth',1.5); hold on;
plot(mo,Otps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])
print('-dpng',[ppath 'Hist_ts_biome_means_zmeso200_schl_sst_noCAN.png'])








