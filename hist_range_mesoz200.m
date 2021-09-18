% CMIP6 output 
% mesoz 200m integrations
% min and max all ESMs and GLMM to get ranges

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_same_orientation.mat');

%%
t(1,1) = nanmin(cz(cz>0));
t(1,2) = nanmax(cz(:));
t(2,1) = nanmin(mz(mz>0));
t(2,2) = nanmax(mz(:));
t(3,1) = nanmin(nz(:));
t(3,2) = nanmax(nz(:));
t(4,1) = nanmin(gz(:));
t(4,2) = nanmax(gz(:));
t(5,1) = nanmin(iz(:));
t(5,2) = nanmax(iz(:));
t(6,1) = nanmin(uz(:));
t(6,2) = nanmax(uz(:));
t(7,1) = nanmin(oz(:));
t(7,2) = nanmax(oz(:));

t(1,3) = quantile(cz(:),0.01);
t(1,4) = quantile(cz(:),0.99);
t(2,3) = quantile(mz(:),0.01);
t(2,4) = quantile(mz(:),0.99);
t(3,3) = quantile(nz(:),0.01);
t(3,4) = quantile(nz(:),0.99);
t(4,3) = quantile(gz(:),0.01);
t(4,4) = quantile(gz(:),0.99);
t(5,3) = quantile(iz(:),0.01);
t(5,4) = quantile(iz(:),0.99);
t(6,3) = quantile(uz(:),0.01);
t(6,4) = quantile(uz(:),0.99);
t(7,3) = quantile(oz(:),0.01);
t(7,4) = quantile(oz(:),0.99);

t(1,5) = quantile(cz(:),0.05);
t(1,6) = quantile(cz(:),0.95);
t(2,5) = quantile(mz(:),0.05);
t(2,6) = quantile(mz(:),0.95);
t(3,5) = quantile(nz(:),0.05);
t(3,6) = quantile(nz(:),0.95);
t(4,5) = quantile(gz(:),0.05);
t(4,6) = quantile(gz(:),0.95);
t(5,5) = quantile(iz(:),0.05);
t(5,6) = quantile(iz(:),0.95);
t(6,5) = quantile(uz(:),0.05);
t(6,6) = quantile(uz(:),0.95);
t(7,5) = quantile(oz(:),0.05);
t(7,6) = quantile(oz(:),0.95);

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
writetable(Trange,[sfile 'range_log10_hist_aclim_200_obsglm.csv'],'WriteRowNames',true);

%% mean and range
apath='/Volumes/MIP/Fish-MIP/CMIP6/';
load([apath 'ISIMIP3b_cellarea_onedeg.mat'])

mr(1,1) = nansum(cz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(2,1) = nansum(mz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(3,1) = nansum(nz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(4,1) = nansum(gz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(5,1) = nansum(iz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(6,1) = nansum(uz(:) .* cell_area(:)) ./ nansum(cell_area(:));
mr(7,1) = nansum(oz(:) .* cell_area(:)) ./ nansum(cell_area(:));

mr(1,2) = quantile(cz(:),0.01);
mr(1,3) = quantile(cz(:),0.99);
mr(2,2) = quantile(mz(:),0.01);
mr(2,3) = quantile(mz(:),0.99);
mr(3,2) = quantile(nz(:),0.01);
mr(3,3) = quantile(nz(:),0.99);
mr(4,2) = quantile(gz(:),0.01);
mr(4,3) = quantile(gz(:),0.99);
mr(5,2) = quantile(iz(:),0.01);
mr(5,3) = quantile(iz(:),0.99);
mr(6,2) = quantile(uz(:),0.01);
mr(6,3) = quantile(uz(:),0.99);
mr(7,2) = quantile(oz(:),0.01);
mr(7,3) = quantile(oz(:),0.99);

mr(:,4) = log10(mr(:,1));
mr(:,5) = log10(mr(:,3)) - log10(mr(:,2));

%%
Mrange = array2table(mr,'RowNames',simtext,'VariableNames',...
    {'areaw_mean','1st','99th','log10areaw_mean','1st_99th'});
writetable(Mrange,[sfile 'mean_range_log10_hist_aclim_200_obsglm.csv'],'WriteRowNames',true);



