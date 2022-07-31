% Calculate ranges of skill metrics

clear all
close all

%%
load('skill_hist_model_obsglm100_climatols.mat')
load('ms_skill_scores_Kendall_hist_model_obsglm100_clim_log10_KK_1e4.mat')

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% remove obsglmm and CAN

bias2 = bias(3:end,:);
stdnorm2 = stdnorm(3:end,:);
tmae2 = tmae(3:end,:);
urmse2 = urmse(3:end,:);

%% model across ESMs
%Pearson
mskill(1,1) = min(Pcc2(:));
mskill(1,2) = mean(Pcc2(:));
mskill(1,3) = max(Pcc2(:));
%Kendall
mskill(2,1) = min(Kcc2(:));
mskill(2,2) = mean(Kcc2(:));
mskill(2,3) = max(Kcc2(:));
%URMSE
mskill(3,1) = min(urmse2(:));
mskill(3,2) = mean(urmse2(:));
mskill(3,3) = max(urmse2(:));
%MAE
mskill(4,1) = min(tmae2(:));
mskill(4,2) = mean(tmae2(:));
mskill(4,3) = max(tmae2(:));
%BIAS
mskill(5,1) = min(abs(bias2(:)));
mskill(5,2) = mean(abs(bias2(:)));
mskill(5,3) = max(abs(bias2(:)));
%NSTD
mskill(6,1) = min(stdnorm2(:));
mskill(6,2) = mean(stdnorm2(:));
mskill(6,3) = max(stdnorm2(:));

%%
Tmskill = array2table(mskill,'RowNames',{'Pearson','Kendall','URMSE','MAE','abs(bias)','normstd'},...
    'VariableNames',{'min','mean','max'});

writetable(Tmskill,[sfile 'ms_mean_skill_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);

%% means by model
%Pearson
eskill(:,1) = mean(Pcc2,2);
%Kendall
eskill(:,2) = mean(Kcc2,2);
%URMSE
eskill(:,3) = mean(urmse2,2);
%MAE
eskill(:,4) = mean(tmae2,2);
%BIAS
eskill(:,5) = mean(abs(bias2),2);
%NSTD
eskill(:,6) = mean(stdnorm2,2);

%%
Teskill = array2table(eskill,'VariableNames',{'Pearson','Kendall','URMSE','MAE','abs(bias)','normstd'},...
    'RowNames',mtex(3:end));

writetable(Teskill,[sfile 'ms_mean_skill_byesm_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);


%%
eskill2=eskill;
%MAE
eskill2(:,4) = mean(10.^tmae2,2);
%BIAS
eskill2(:,5) = mean(abs(10.^bias2),2);










