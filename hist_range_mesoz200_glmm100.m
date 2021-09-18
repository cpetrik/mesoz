% CMIP6 output 
% mesoz 200m integrations
% min and max all ESMs and GLMM to get ranges

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat');

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
czmo_all = czmo_all * 12.01 * 1e3;
mzmo_all = mzmo_all * 12.01 * 1e3;
nzmo_all = nzmo_all * 12.01 * 1e3;
gzmo_all = gzmo_all * 12.01 * 1e3;
izmo_all = izmo_all * 12.01 * 1e3;
uzmo_all = uzmo_all * 12.01 * 1e3;

%%
t(1,1) = nanmin(czmo_all(czmo_all>0));
t(1,2) = nanmax(czmo_all(:));
t(2,1) = nanmin(mzmo_all(mzmo_all>0));
t(2,2) = nanmax(mzmo_all(:));
t(3,1) = nanmin(nzmo_all(:));
t(3,2) = nanmax(nzmo_all(:));
t(4,1) = nanmin(gzmo_all(:));
t(4,2) = nanmax(gzmo_all(:));
t(5,1) = nanmin(izmo_all(:));
t(5,2) = nanmax(izmo_all(:));
t(6,1) = nanmin(uzmo_all(:));
t(6,2) = nanmax(uzmo_all(:));
t(7,1) = nanmin(ozmo_all(:));
t(7,2) = nanmax(ozmo_all(:));

t(1,3) = quantile(czmo_all(:),0.01);
t(1,4) = quantile(czmo_all(:),0.99);
t(2,3) = quantile(mzmo_all(:),0.01);
t(2,4) = quantile(mzmo_all(:),0.99);
t(3,3) = quantile(nzmo_all(:),0.01);
t(3,4) = quantile(nzmo_all(:),0.99);
t(4,3) = quantile(gzmo_all(:),0.01);
t(4,4) = quantile(gzmo_all(:),0.99);
t(5,3) = quantile(izmo_all(:),0.01);
t(5,4) = quantile(izmo_all(:),0.99);
t(6,3) = quantile(uzmo_all(:),0.01);
t(6,4) = quantile(uzmo_all(:),0.99);
t(7,3) = quantile(ozmo_all(:),0.01);
t(7,4) = quantile(ozmo_all(:),0.99);

t(1,5) = quantile(czmo_all(:),0.05);
t(1,6) = quantile(czmo_all(:),0.95);
t(2,5) = quantile(mzmo_all(:),0.05);
t(2,6) = quantile(mzmo_all(:),0.95);
t(3,5) = quantile(nzmo_all(:),0.05);
t(3,6) = quantile(nzmo_all(:),0.95);
t(4,5) = quantile(gzmo_all(:),0.05);
t(4,6) = quantile(gzmo_all(:),0.95);
t(5,5) = quantile(izmo_all(:),0.05);
t(5,6) = quantile(izmo_all(:),0.95);
t(6,5) = quantile(uzmo_all(:),0.05);
t(6,6) = quantile(uzmo_all(:),0.95);
t(7,5) = quantile(ozmo_all(:),0.05);
t(7,6) = quantile(ozmo_all(:),0.95);

%%
t10 = log10(t);
t10(:,7) = t10(:,2) - t10(:,1);
t10(:,8) = t10(:,4) - t10(:,3);
t10(:,9) = t10(:,6) - t10(:,5);

simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK','obs'};
sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

range = t10(:,7:9);
Trange = array2table(range,'RowNames',simtext,'VariableNames',...
    {'full','1st_99th','5th_95th'});
writetable(Trange,[sfile 'range_log10_hist_aclim_200_obsglm100.csv'],'WriteRowNames',true);

%% mean and range
apath='/Volumes/MIP/Fish-MIP/CMIP6/';
load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

mr(1,1) = nansum(czmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(2,1) = nansum(mzmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(3,1) = nansum(nzmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(4,1) = nansum(gzmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(5,1) = nansum(izmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(6,1) = nansum(uzmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(7,1) = nansum(ozmo_all(:) .* cell_area(:)) ./ nansum(cell_area(:));

mr(1,2) = quantile(czmo_all(:),0.01);
mr(1,3) = quantile(czmo_all(:),0.99);
mr(2,2) = quantile(mzmo_all(:),0.01);
mr(2,3) = quantile(mzmo_all(:),0.99);
mr(3,2) = quantile(nzmo_all(:),0.01);
mr(3,3) = quantile(nzmo_all(:),0.99);
mr(4,2) = quantile(gzmo_all(:),0.01);
mr(4,3) = quantile(gzmo_all(:),0.99);
mr(5,2) = quantile(izmo_all(:),0.01);
mr(5,3) = quantile(izmo_all(:),0.99);
mr(6,2) = quantile(uzmo_all(:),0.01);
mr(6,3) = quantile(uzmo_all(:),0.99);
mr(7,2) = quantile(ozmo_all(:),0.01);
mr(7,3) = quantile(ozmo_all(:),0.99);

mr(:,4) = log10(mr(:,1));
mr(:,5) = log10(mr(:,3)) - log10(mr(:,2));

%%
Mrange = array2table(mr,'RowNames',simtext,'VariableNames',...
    {'areaw_mean','1st','99th','log10areaw_mean','1st_99th'});
writetable(Mrange,[sfile 'mean_range_log10_hist_aclim_200_obsglm100.csv'],'WriteRowNames',true);


