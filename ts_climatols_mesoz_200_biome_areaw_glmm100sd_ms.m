% CMIP6 output 
% 200m integrations
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('Hist_ts_chl_clim_biomes_areaw.mat')
load('Hist_ts_sst_clim_biomes_areaw.mat')
load('Hist_ts_mesoz_clim_biomes_areaw.mat')
load('glm100_output_by_biom_v2.mat')

%% Plots 

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

%% Std devs as shadded region
%create continuous x value array for plotting
mo = 1:12;

X=[mo fliplr(mo)]; 
%create y values for out and then back
%+/- 2 stdev
Sz=[OGmzmean+2*OGmzsd; flipud(OGmzmean-2*OGmzsd)]; 

%% Try using tiles
close all

figure(1)
tiledlayout(4,3, 'TileSpacing', 'compact')
% MESOZ ------------------------------
%subplot(4,3,1)
nexttile
plot(mo,Cz,'LineWidth',1.5); hold on;
plot(mo,Mz,'LineWidth',1.5); hold on;
plot(mo,Nz,'LineWidth',1.5); hold on;
plot(mo,Gz,'LineWidth',1.5); hold on;
plot(mo,Iz,'LineWidth',1.5); hold on;
plot(mo,Uz,'LineWidth',1.5); hold on;
%plot(mo,Oz,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,1),'LineWidth',1.5); hold on;
fill(X,Sz(:,1),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('mesozoo (mgC m^-^2)')
ylabel('Global')
xlim([1 12])
ylim([0 1500])
text(1,1970,'A','FontWeight','Bold','FontSize',14)

% CHL --------------------------------
nexttile %subplot(4,3,2)
plot(mo,Cc,'LineWidth',1.5); hold on;
plot(mo,Mc,'LineWidth',1.5); hold on;
plot(mo,Nc,'LineWidth',1.5); hold on;
plot(mo,Gc,'LineWidth',1.5); hold on;
plot(mo,Ic,'LineWidth',1.5); hold on;
plot(mo,Uc,'LineWidth',1.5); hold on;
plot(mo,Oc,'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
ylabel('Global')
xlim([1 12])
text(1,0.525,'B','FontWeight','Bold','FontSize',14)

% SST --------------------------
nexttile %subplot(4,3,3)
plot(mo,Ct,'LineWidth',1.5); hold on;
plot(mo,Mt,'LineWidth',1.5); hold on;
plot(mo,Nt,'LineWidth',1.5); hold on;
plot(mo,Gt,'LineWidth',1.5); hold on;
plot(mo,It,'LineWidth',1.5); hold on;
plot(mo,Ut,'LineWidth',1.5); hold on;
plot(mo,Ot,'LineWidth',1.5); hold on;
title('SST (^oC)')
ylabel('Global')
xlim([1 12])
text(1,16.5,'C','FontWeight','Bold','FontSize',14)

nexttile %(4,3,4)
plot(mo,Czlc,'LineWidth',1.5); hold on;
plot(mo,Mzlc,'LineWidth',1.5); hold on;
plot(mo,Nzlc,'LineWidth',1.5); hold on;
plot(mo,Gzlc,'LineWidth',1.5); hold on;
plot(mo,Izlc,'LineWidth',1.5); hold on;
plot(mo,Uzlc,'LineWidth',1.5); hold on;
%plot(mo,Ozlc,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,2),'LineWidth',1.5); hold on;
fill(X,Sz(:,2),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,5)
plot(mo,Cclc,'LineWidth',1.5); hold on;
plot(mo,Mclc,'LineWidth',1.5); hold on;
plot(mo,Nclc,'LineWidth',1.5); hold on;
plot(mo,Gclc,'LineWidth',1.5); hold on;
plot(mo,Iclc,'LineWidth',1.5); hold on;
plot(mo,Uclc,'LineWidth',1.5); hold on;
plot(mo,Oclc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,6)
plot(mo,Ctlc,'LineWidth',1.5); hold on;
plot(mo,Mtlc,'LineWidth',1.5); hold on;
plot(mo,Ntlc,'LineWidth',1.5); hold on;
plot(mo,Gtlc,'LineWidth',1.5); hold on;
plot(mo,Itlc,'LineWidth',1.5); hold on;
plot(mo,Utlc,'LineWidth',1.5); hold on;
plot(mo,Otlc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %(4,3,7)
plot(mo,Czss,'LineWidth',1.5); hold on;
plot(mo,Mzss,'LineWidth',1.5); hold on;
plot(mo,Nzss,'LineWidth',1.5); hold on;
plot(mo,Gzss,'LineWidth',1.5); hold on;
plot(mo,Izss,'LineWidth',1.5); hold on;
plot(mo,Uzss,'LineWidth',1.5); hold on;
%plot(mo,Ozss,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,3),'LineWidth',1.5); hold on;
fill(X,Sz(:,3),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,8)
plot(mo,Ccss,'LineWidth',1.5); hold on;
plot(mo,Mcss,'LineWidth',1.5); hold on;
plot(mo,Ncss,'LineWidth',1.5); hold on;
plot(mo,Gcss,'LineWidth',1.5); hold on;
plot(mo,Icss,'LineWidth',1.5); hold on;
plot(mo,Ucss,'LineWidth',1.5); hold on;
plot(mo,Ocss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,9)
plot(mo,Ctss,'LineWidth',1.5); hold on;
plot(mo,Mtss,'LineWidth',1.5); hold on;
plot(mo,Ntss,'LineWidth',1.5); hold on;
plot(mo,Gtss,'LineWidth',1.5); hold on;
plot(mo,Itss,'LineWidth',1.5); hold on;
plot(mo,Utss,'LineWidth',1.5); hold on;
plot(mo,Otss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %(4,3,10)
plot(mo,Czps,'LineWidth',1.5); hold on;
plot(mo,Mzps,'LineWidth',1.5); hold on;
plot(mo,Nzps,'LineWidth',1.5); hold on;
plot(mo,Gzps,'LineWidth',1.5); hold on;
plot(mo,Izps,'LineWidth',1.5); hold on;
plot(mo,Uzps,'LineWidth',1.5); hold on;
%plot(mo,Ozps,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,4),'LineWidth',1.5); hold on;
fill(X,Sz(:,4),'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('HCPS')
xlim([1 12])

nexttile %subplot(4,3,11)
plot(mo,Ccps,'LineWidth',1.5); hold on;
plot(mo,Mcps,'LineWidth',1.5); hold on;
plot(mo,Ncps,'LineWidth',1.5); hold on;
plot(mo,Gcps,'LineWidth',1.5); hold on;
plot(mo,Icps,'LineWidth',1.5); hold on;
plot(mo,Ucps,'LineWidth',1.5); hold on;
plot(mo,Ocps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlabel('NH month')
xlim([1 12])

nexttile %subplot(4,3,12)
plot(mo,Ctps,'LineWidth',1.5); hold on;
plot(mo,Mtps,'LineWidth',1.5); hold on;
plot(mo,Ntps,'LineWidth',1.5); hold on;
plot(mo,Gtps,'LineWidth',1.5); hold on;
plot(mo,Itps,'LineWidth',1.5); hold on;
plot(mo,Utps,'LineWidth',1.5); hold on;
plot(mo,Otps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])

lg = legend(nexttile(3),{'CAN','CMCC','CNRM','GFDL','IPSL','UKESM','obsGLMM/sat'});
lg.Location = 'eastoutside';
print('-dpng',[ppath 'Hist_clim_biome_areaw_means_zmeso200_schl_sst_glmm100SDs.png'])

%% No obs chl, sst

figure(2)
tiledlayout(4,3, 'TileSpacing', 'compact')
% MESOZ ------------------------------
%subplot(4,3,1)
nexttile
plot(mo,Cz,'LineWidth',1.5); hold on;
plot(mo,Mz,'LineWidth',1.5); hold on;
plot(mo,Nz,'LineWidth',1.5); hold on;
plot(mo,Gz,'LineWidth',1.5); hold on;
plot(mo,Iz,'LineWidth',1.5); hold on;
plot(mo,Uz,'LineWidth',1.5); hold on;
plot(mo,Oz,'LineWidth',1.5); hold on;
title('mesozoo (mgC m^-^2)')
ylabel('Global')
xlim([1 12])

% CHL --------------------------------
nexttile %subplot(4,3,2)
plot(mo,Cc,'LineWidth',1.5); hold on;
plot(mo,Mc,'LineWidth',1.5); hold on;
plot(mo,Nc,'LineWidth',1.5); hold on;
plot(mo,Gc,'LineWidth',1.5); hold on;
plot(mo,Ic,'LineWidth',1.5); hold on;
plot(mo,Uc,'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
ylabel('Global')
xlim([1 12])

% SST --------------------------
nexttile %subplot(4,3,3)
plot(mo,Ct,'LineWidth',1.5); hold on;
plot(mo,Mt,'LineWidth',1.5); hold on;
plot(mo,Nt,'LineWidth',1.5); hold on;
plot(mo,Gt,'LineWidth',1.5); hold on;
plot(mo,It,'LineWidth',1.5); hold on;
plot(mo,Ut,'LineWidth',1.5); hold on;
title('SST (^oC)')
ylabel('Global')
xlim([1 12])

nexttile %(4,3,4)
plot(mo,Czlc,'LineWidth',1.5); hold on;
plot(mo,Mzlc,'LineWidth',1.5); hold on;
plot(mo,Nzlc,'LineWidth',1.5); hold on;
plot(mo,Gzlc,'LineWidth',1.5); hold on;
plot(mo,Izlc,'LineWidth',1.5); hold on;
plot(mo,Uzlc,'LineWidth',1.5); hold on;
plot(mo,Ozlc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,5)
plot(mo,Cclc,'LineWidth',1.5); hold on;
plot(mo,Mclc,'LineWidth',1.5); hold on;
plot(mo,Nclc,'LineWidth',1.5); hold on;
plot(mo,Gclc,'LineWidth',1.5); hold on;
plot(mo,Iclc,'LineWidth',1.5); hold on;
plot(mo,Uclc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,6)
plot(mo,Ctlc,'LineWidth',1.5); hold on;
plot(mo,Mtlc,'LineWidth',1.5); hold on;
plot(mo,Ntlc,'LineWidth',1.5); hold on;
plot(mo,Gtlc,'LineWidth',1.5); hold on;
plot(mo,Itlc,'LineWidth',1.5); hold on;
plot(mo,Utlc,'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %(4,3,7)
plot(mo,Czss,'LineWidth',1.5); hold on;
plot(mo,Mzss,'LineWidth',1.5); hold on;
plot(mo,Nzss,'LineWidth',1.5); hold on;
plot(mo,Gzss,'LineWidth',1.5); hold on;
plot(mo,Izss,'LineWidth',1.5); hold on;
plot(mo,Uzss,'LineWidth',1.5); hold on;
plot(mo,Ozss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,8)
plot(mo,Ccss,'LineWidth',1.5); hold on;
plot(mo,Mcss,'LineWidth',1.5); hold on;
plot(mo,Ncss,'LineWidth',1.5); hold on;
plot(mo,Gcss,'LineWidth',1.5); hold on;
plot(mo,Icss,'LineWidth',1.5); hold on;
plot(mo,Ucss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,9)
plot(mo,Ctss,'LineWidth',1.5); hold on;
plot(mo,Mtss,'LineWidth',1.5); hold on;
plot(mo,Ntss,'LineWidth',1.5); hold on;
plot(mo,Gtss,'LineWidth',1.5); hold on;
plot(mo,Itss,'LineWidth',1.5); hold on;
plot(mo,Utss,'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %(4,3,10)
plot(mo,Czps,'LineWidth',1.5); hold on;
plot(mo,Mzps,'LineWidth',1.5); hold on;
plot(mo,Nzps,'LineWidth',1.5); hold on;
plot(mo,Gzps,'LineWidth',1.5); hold on;
plot(mo,Izps,'LineWidth',1.5); hold on;
plot(mo,Uzps,'LineWidth',1.5); hold on;
plot(mo,Ozps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])

nexttile %subplot(4,3,11)
plot(mo,Ccps,'LineWidth',1.5); hold on;
plot(mo,Mcps,'LineWidth',1.5); hold on;
plot(mo,Ncps,'LineWidth',1.5); hold on;
plot(mo,Gcps,'LineWidth',1.5); hold on;
plot(mo,Icps,'LineWidth',1.5); hold on;
plot(mo,Ucps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlabel('NH month')
xlim([1 12])

nexttile %subplot(4,3,12)
plot(mo,Ctps,'LineWidth',1.5); hold on;
plot(mo,Mtps,'LineWidth',1.5); hold on;
plot(mo,Ntps,'LineWidth',1.5); hold on;
plot(mo,Gtps,'LineWidth',1.5); hold on;
plot(mo,Itps,'LineWidth',1.5); hold on;
plot(mo,Utps,'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])

lg = legend(nexttile(1),{'CAN','CMCC','CNRM','GFDL','IPSL','UKESM','obsGLMM'});
lg.Location = 'westoutside';
print('-dpng',[ppath 'Hist_clim_biome_areaw_means_zmeso200_schl_sst_glmm100SDs_noOBS.png'])







