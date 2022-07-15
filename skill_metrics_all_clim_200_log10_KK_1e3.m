% Calculate different skill metrics for each ESM
% log10 transform biomass
% uses SH shift
% use Kelly Kearney's skill cal fn based on Stowe et a.

clear all
close all

%%
load('skill_hist_model_obsglm100_climatols.mat')

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% Standardization
% use quantiles to determine what value to change zeros?
quantile(comb(:,3:9),[0.01 0.05])
% limit of observing capabilities 1 mgC/m2

%% log10
Acomb = log10(comb(:,3:9)+1e-4);
Dcomb = log10(dcomb(:,3:9)+1e-4);
Jcomb = log10(jcomb(:,3:9)+1e-4);
Mcomb = log10(mcomb(:,3:9)+1e-4);
Scomb = log10(scomb(:,3:9)+1e-4);

%% Skill metric = weighted sum of squares for multivariate
metrics={'std','r','RMSD','CRMSD','bias','stdnorm','rmsdnorm','crmsdnorm',...
    'AAE','RI','MEF'};
metrics=metrics';

%   S:  1 x 1 structure with the following fields; all are 1 x m+1 arrays,
%       where the first element corresponds to the reference data and 2:end
%       correspond to the models.
%       std:        standard deviation, normalized to n.
%       cor:        correlation coefficient
%       rmsd:       root mean squared difference
%       crmsd:      centered pattern (i.e. unbiased) root mean squared
%                   difference
%       bias:       bias, i.e. average error
%       stdnorm:    normalized standard deviation
%       rmsdnorm:   root mean squared difference, normalized to standard
%                   deviation of reference data
%       crmsdnorm:  centered root mean squared difference, normalized to
%                   standard deviation of reference data
%       aae:        average absolute error
%       ri:         reliability index (factor by which model differs from
%                   reference... note that this metric falls apart if any
%                   of the modeled values are 0).
%       mef:        modeling efficiency (skill relative to average of
%                   observations, 1 = perfect, 0 = same as averaging obs,
%                   <1 = worse than just averaging observations)

Amod_all = Acomb;
o=(Amod_all(:,7));
m=(Amod_all(:,1:6));
A = skillstats2(o, m);

Dmod_all = Dcomb;
o=(Dmod_all(:,7));
m=(Dmod_all(:,1:6));
D = skillstats2(o, m);

Mmod_all = Mcomb;
o=(Mmod_all(:,7));
m=(Mmod_all(:,1:6));
M = skillstats2(o, m);

Jmod_all = Jcomb;
o=(Jmod_all(:,7));
m=(Jmod_all(:,1:6));
J = skillstats2(o, m);

Smod_all = Scomb;
o=(Smod_all(:,7));
m=(Smod_all(:,1:6));
S = skillstats2(o, m);

%%
skill=NaN*ones(11,7,5);
skill(1,:,1) = A.std;
skill(2,:,1) = A.cor;
skill(3,:,1) = A.rmsd;
skill(4,:,1) = A.crmsd;
skill(5,:,1) = A.bias;
skill(6,:,1) = A.stdnorm;
skill(7,:,1) = A.rmsdnorm;
skill(8,:,1) = A.crmsdnorm;
skill(9,:,1) = A.aae;
skill(10,:,1) = A.ri;
skill(11,:,1) = A.mef;

skill(1,:,2) = D.std;
skill(2,:,2) = D.cor;
skill(3,:,2) = D.rmsd;
skill(4,:,2) = D.crmsd;
skill(5,:,2) = D.bias;
skill(6,:,2) = D.stdnorm;
skill(7,:,2) = D.rmsdnorm;
skill(8,:,2) = D.crmsdnorm;
skill(9,:,2) = D.aae;
skill(10,:,2) = D.ri;
skill(11,:,2) = D.mef;

skill(1,:,3) = M.std;
skill(2,:,3) = M.cor;
skill(3,:,3) = M.rmsd;
skill(4,:,3) = M.crmsd;
skill(5,:,3) = M.bias;
skill(6,:,3) = M.stdnorm;
skill(7,:,3) = M.rmsdnorm;
skill(8,:,3) = M.crmsdnorm;
skill(9,:,3) = M.aae;
skill(10,:,3) = M.ri;
skill(11,:,3) = M.mef;

skill(1,:,4) = J.std;
skill(2,:,4) = J.cor;
skill(3,:,4) = J.rmsd;
skill(4,:,4) = J.crmsd;
skill(5,:,4) = J.bias;
skill(6,:,4) = J.stdnorm;
skill(7,:,4) = J.rmsdnorm;
skill(8,:,4) = J.crmsdnorm;
skill(9,:,4) = J.aae;
skill(10,:,4) = J.ri;
skill(11,:,4) = J.mef;

skill(1,:,5) = S.std;
skill(2,:,5) = S.cor;
skill(3,:,5) = S.rmsd;
skill(4,:,5) = S.crmsd;
skill(5,:,5) = S.bias;
skill(6,:,5) = S.stdnorm;
skill(7,:,5) = S.rmsdnorm;
skill(8,:,5) = S.crmsdnorm;
skill(9,:,5) = S.aae;
skill(10,:,5) = S.ri;
skill(11,:,5) = S.mef;

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey
% set(groot,'defaultAxesColorOrder',cm);

simtext = {'obsGLMM','CAN','CMCC','CNRM','GFDL','IPSL','UK'};
simtex = {'CA','CM','CN','GF','IP','UK'};
ctext = {'All','Winter','Spring','Summer','Fall'};

save([sfile 'skill_scores_hist_model_obsglm100_clim_log10_KK_1e4.mat'],'skill',...
    'simtext','metrics','ctext')

%% Tables
cor = squeeze(skill(2,:,:));
Tcorr = array2table(cor,'RowNames',simtext,'VariableNames',ctext);

rmsd = squeeze(skill(3,:,:));
Trmse = array2table(rmsd,'RowNames',simtext,'VariableNames',ctext);

stdnorm = squeeze(skill(6,:,:));
Tnstd = array2table(stdnorm,'RowNames',simtext,'VariableNames',ctext);

urmse = squeeze(skill(4,:,:));
Turmse = array2table(urmse,'RowNames',simtext,'VariableNames',ctext);

tmae = squeeze(skill(9,:,:));
Tmae = array2table(tmae,'RowNames',simtext,'VariableNames',ctext);

bias = squeeze(skill(5,:,:));
Tbias = array2table(bias,'RowNames',simtext,'VariableNames',ctext);

writetable(Tcorr,[sfile 'corr_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);
writetable(Trmse,[sfile 'rmse_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);
writetable(Tnstd,[sfile 'nstd_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);
writetable(Turmse,[sfile 'urmse_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);
writetable(Tmae,[sfile 'mae_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);
writetable(Tbias,[sfile 'bias_hist_clims_200_obsglm100_log10_KK_1e4.csv'],'WriteRowNames',true);

%% make data vecs into tables to also calc stats in R

TAcomb = array2table(Acomb,'VariableNames',simtext);
TDcomb = array2table(Dcomb,'VariableNames',simtext);
TJcomb = array2table(Jcomb,'VariableNames',simtext);
TMcomb = array2table(Mcomb,'VariableNames',simtext);
TScomb = array2table(Scomb,'VariableNames',simtext);

writetable(TAcomb,[sfile 'climatol_All_hist_clims_200_obsglm100_log10_1e4.csv'],'WriteRowNames',false);
writetable(TDcomb,[sfile 'climatol_DJF_hist_clims_200_obsglm100_log10_1e4.csv'],'WriteRowNames',false);
writetable(TJcomb,[sfile 'climatol_JJA_hist_clims_200_obsglm100_log10_1e4.csv'],'WriteRowNames',false);
writetable(TMcomb,[sfile 'climatol_MAM_hist_clims_200_obsglm100_log10_1e4.csv'],'WriteRowNames',false);
writetable(TScomb,[sfile 'climatol_SON_hist_clims_200_obsglm100_log10_1e4.csv'],'WriteRowNames',false);

