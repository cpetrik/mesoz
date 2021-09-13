% CMIP6 output 
% 200m integrations
% Area-weighted

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% not area-weighted
load('Hist_ts_obs_clim_biomes.mat')
Oc2 = Oc;
Oclc2 = Oclc;
Ocss2 = Ocss;
Ocps2 = Ocps;
Ot2 = Ot;
Otlc2 = Otlc;
Otss2 = Otss;
Otps2 = Otps;
Oz2 = Oz;
Ozlc2 = Ozlc;
Ozss2 = Ozss;
Ozps2 = Ozps;

%%
load('Hist_ts_chl_clim_biomes_areaw.mat')
load('Hist_ts_sst_clim_biomes_areaw.mat')
load('Hist_ts_mesoz_clim_biomes_areaw.mat')
load('glm_output_by_biom_v2.mat')

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
plot(mo,Oz2,'LineWidth',1.5); hold on;
plot(mo,Oz,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,1),'LineWidth',1.5); hold on;
fill(X,Sz(:,1),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('mesozoo (mgC m^-^2)')
ylabel('Global')
xlim([1 12])

% CHL --------------------------------
nexttile %subplot(4,3,2)
plot(mo,Oc2,'LineWidth',1.5); hold on;
plot(mo,Oc,'LineWidth',1.5); hold on;
plot(mo,OGchl(:,1),'LineWidth',1.5); hold on;
title('chl (mg m^-^3)')
ylabel('Global')
xlim([1 12])

% SST --------------------------
nexttile %subplot(4,3,3)
plot(mo,Ot2,'LineWidth',1.5); hold on;
plot(mo,Ot,'LineWidth',1.5); hold on;
plot(mo,OGsst(:,1),'LineWidth',1.5); hold on;
title('SST (^oC)')
ylabel('Global')
xlim([1 12])

nexttile %(4,3,4)
plot(mo,Ozlc2,'LineWidth',1.5); hold on;
plot(mo,Ozlc,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,2),'LineWidth',1.5); hold on;
fill(X,Sz(:,2),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,5)
plot(mo,Oclc2,'LineWidth',1.5); hold on;
plot(mo,Oclc,'--','LineWidth',1.5); hold on;
plot(mo,OGchl(:,2),'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %subplot(4,3,6)
plot(mo,Otlc2,'LineWidth',1.5); hold on;
plot(mo,Otlc,'LineWidth',1.5); hold on;
plot(mo,OGsst(:,2),'LineWidth',1.5); hold on;
ylabel('LC')
xlim([1 12])

nexttile %(4,3,7)
plot(mo,Ozss2,'LineWidth',1.5); hold on;
plot(mo,Ozss,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,3),'LineWidth',1.5); hold on;
fill(X,Sz(:,3),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,8)
plot(mo,Ocss2,'LineWidth',1.5); hold on;
plot(mo,Ocss,'LineWidth',1.5); hold on;
plot(mo,OGchl(:,3),'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %subplot(4,3,9)
plot(mo,Otss2,'LineWidth',1.5); hold on;
plot(mo,Otss,'LineWidth',1.5); hold on;
plot(mo,OGsst(:,3),'LineWidth',1.5); hold on;
ylabel('HCSS')
xlim([1 12])

nexttile %(4,3,10)
plot(mo,Ozps2,'LineWidth',1.5); hold on;
plot(mo,Ozps,'LineWidth',1.5); hold on;
plot(mo,OGmzmean(:,4),'LineWidth',1.5); hold on;
fill(X,Sz(:,4),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylabel('HCPS')
xlim([1 12])

nexttile %subplot(4,3,11)
plot(mo,Ocps2,'LineWidth',1.5); hold on;
plot(mo,Ocps,'LineWidth',1.5); hold on;
plot(mo,OGchl(:,4),'LineWidth',1.5); hold on;
ylabel('HCPS')
xlabel('NH month')
xlim([1 12])

nexttile %subplot(4,3,12)
plot(mo,Otps2,'LineWidth',1.5); hold on;
plot(mo,Otps,'LineWidth',1.5); hold on;
plot(mo,OGsst(:,4),'LineWidth',1.5); hold on;
ylabel('HCPS')
xlim([1 12])

lg = legend(nexttile(3),{'ColleenNW','ColleenW','RyanW'});
lg.Location = 'eastoutside';
print('-dpng',[ppath 'Comp_Hist_clim_biome_areaw_means_zmeso200_schl_sst_SDs.png'])






