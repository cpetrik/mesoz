% Calculate different skill metrics for each ESM
% log10 transform biomass
% uses SH shift

clear all
close all

%%
load('skill_hist_model_obsglm100_climatols.mat')

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% No Standardization
% log10
Acomb = log10(comb(:,3:9)+1);
Dcomb = log10(dcomb(:,3:9)+1);
Jcomb = log10(jcomb(:,3:9)+1);
Mcomb = log10(mcomb(:,3:9)+1);
Scomb = log10(scomb(:,3:9)+1);

%% Skill metric = weighted sum of squares for multivariate
metrics={'ccoef','crmsd','sdev'};
metrics=metrics';

skill=NaN*ones(3,6,5);
for j=1:5
    for k=1:6
        % Pick
        if j==1
            mod_all = Acomb;
        elseif j==2
            mod_all = Dcomb;
        elseif j==3
            mod_all = Mcomb;
        elseif j==4
            mod_all = Jcomb;
        else
            mod_all = Scomb;
        end
        
        o=(mod_all(:,7));
        p=(mod_all(:,k));
        
        % Only compare model & obs with non-NaNs
        p(isnan(o))=[];
        o(isnan(o))=[];
        o(isnan(p))=[];
        p(isnan(p))=[];
        
        [stats] = taylor_statistics(p,o);
        skill(1,k,j) = stats.ccoef(2);
        skill(2,k,j) = stats.crmsd(2);
        skill(3,k,j) = stats.sdev(2);
        
    end
end

simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK'};
ctext = {'All','Winter','Spring','Summer','Fall'};

save('skill_scores_hist_model_obsglm100_clim_log10_taylor.mat','skill','simtext','metrics','ctext')

