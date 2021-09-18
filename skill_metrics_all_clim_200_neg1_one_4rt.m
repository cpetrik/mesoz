% Calculate different skill metrics for each ESM
% 4th-rt transform biomass
% uses SH shift

clear all
close all

%%
load('skill_hist_model_obsglm100_climatols.mat')

sfile = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/data_stats_zmeso/';

%% Standardization
% 1st take sqrt
Acomb = (comb(:,3:9)).^(1/4);
Dcomb = (dcomb(:,3:9)).^(1/4);
Jcomb = (jcomb(:,3:9)).^(1/4);
Mcomb = (mcomb(:,3:9)).^(1/4);
Scomb = (scomb(:,3:9)).^(1/4);

% then stdize to -1 to 1
Amean = nanmean(Acomb);
Astd  = nanstd(Acomb);
Acomb = (Acomb - Amean) ./ Astd;

Dmean = nanmean(Dcomb);
Dstd  = nanstd(Dcomb);
Dcomb = (Dcomb - Dmean) ./ Dstd;

Jmean = nanmean(Jcomb);
Jstd  = nanstd(Jcomb);
Jcomb = (Jcomb - Jmean) ./ Jstd;

Mmean = nanmean(Mcomb);
Mstd  = nanstd(Mcomb);
Mcomb = (Mcomb - Mmean) ./ Mstd;

Smean = nanmean(Scomb);
Sstd  = nanstd(Scomb);
Scomb = (Scomb - Smean) ./ Sstd;

%% Skill metric = weighted sum of squares for multivariate
metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

skill=NaN*ones(14,6,5);
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
        
        n = length(o);
        
        omean=repmat(nanmean(o),n,1);
        pmean=repmat(nanmean(p),n,1);
        osig=nanstd(o);
        psig=nanstd(p);
        
        % corr coeff
        num=nansum((o-omean).*(p-pmean));
        d1=nansum((o-omean).^2);
        d2=nansum((p-pmean).^2);
        den=sqrt(d1*d2);
        skill(1,k,j) = num/den;
        
        % root mean square error
        num=nansum((p-o).^2);
        skill(2,k,j) = sqrt(num/n);
        
        % reliability index  %%What about predicted zero values?
        %q1=nansum((log(o./p)).^2);
        q1=nansum((o-p).^2);
        skill(3,k,j) = exp(sqrt(q1/n));
        
        % average error
        skill(4,k,j) = nansum(p-o) / n;
        
        % average absolute error
        skill(5,k,j) = nansum(abs(p-o)) / n;
        
        % modeling efficiency
        num1=nansum((o-omean).^2);
        num2=nansum((p-o).^2);
        skill(6,k,j) = (num1-num2)/num1;
        
        % Taylor R
        num=nansum((o-omean).*(p-pmean));
        skill(7,k,j) = num/(n*osig*psig);
        
        % Taylor normalized std
        skill(8,k,j) = psig/osig;
        
        % unbiased root mean square difference
        % sign tells difference between model and obs std
        q1=nansum(((p-pmean)-(o-omean)).^2);
        if (psig>=osig)
            skill(9,k,j) = sqrt(q1/n);
        else
            skill(9,k,j) = -1*sqrt(q1/n);
        end
        
        % total root mean square difference
        q1=nansum((o-p).^2);
        skill(10,k,j) = sqrt(q1/n);
        
        % normalized bias
        skill(11,k,j) = (pmean(1,1)-omean(1,1));%./osig;
        
        % Joliff et al 2009 S1
        R=skill(7,k,j);
        s=skill(8,k,j);
        num=2*(1+R);
        den=(s + (1/s)).^2;
        skill(12,k,j) = 1 - (num/den);
        
        % Joliff et al 2009 S2
        num=(1+R).^4;
        den=4*(s + (1/s)).^2;
        skill(13,k,j) = 1 - (num/den);
        
        % Joliff et al 2009 S3
        q1=exp(-((s-1).^2) / 0.18);
        q2=(1+R)/2;
        skill(14,k,j) = 1 - (q1*q2);
    end
end

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey
% set(groot,'defaultAxesColorOrder',cm);

simtext = {'CAN','CMCC','CNRM','GFDL','IPSL','UK'};
simtex = {'CA','CM','CN','GF','IP','UK'};
ctext = {'All','Winter','Spring','Summer','Fall'};

save('skill_scores_hist_model_obsglm100_clim_neg11_4rt.mat','skill','simtext','metrics','ctext')

% Bar graphs
% figure(1)
% subplot(3,5,1)
% bar(skill(1,:,1),'k')
% title('Annual Corr coeff')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,2)
% bar(skill(1,:,2),'k')
% title('Winter')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,3)
% bar(skill(1,:,3),'k')
% title('Spring')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,4)
% bar(skill(1,:,4),'k')
% title('Summer')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,5)
% bar(skill(1,:,5),'k')
% title('Fall')
% set(gca,'XTickLabel',simtex)
% 
% 
% subplot(3,5,6)
% bar(skill(2,:,4),'k')
% %title('Root mean square error')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,7)
% bar(skill(2,:,1),'k')
% %title('Root mean square error')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,8)
% bar(skill(2,:,2),'k')
% title('Root mean square error')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,9)
% bar(skill(2,:,3),'k')
% %title('Root mean square error')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,10)
% bar(skill(2,:,5),'k')
% %title('Root mean square error')
% set(gca,'XTickLabel',simtex)
% 
% 
% subplot(3,5,11)
% bar(skill(11,:,1),'k')
% %title('Bias')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,12)
% bar(skill(11,:,2),'k')
% %title('Bias')
% set(gca,'XTickLabel',simtex)
% 
% subplot(3,5,13)
% bar(skill(11,:,3),'k')
% title('Bias')
% set(gca,'XTickLabel',simtex)
% %ylim([-15 5])
% 
% subplot(3,5,14)
% bar(skill(11,:,4),'k')
% %title('Bias')
% set(gca,'XTickLabel',simtex)
% %ylim([-15 5])
% 
% subplot(3,5,15)
% bar(skill(11,:,5),'k')
% %title('Bias')
% set(gca,'XTickLabel',simtex)
% %ylim([-15 5])
% print('-dpng',[figp 'corr_rmse_bias_clims_200_obsglm_neg11_4rt_subplot.png'])

%%

% figure(6)
% subplot(3,2,1)
% bar(squeeze(skill(1,:,:))')
% title('Correlation coefficient')
% set(gca,'XTickLabel',ctext)
% 
% subplot(3,2,3)
% bar(squeeze(skill(2,:,:))')
% title('Root mean square error')
% set(gca,'XTickLabel',ctext)
% 
% subplot(3,2,5)
% bar(squeeze(skill(11,:,:))')
% title('Bias')
% set(gca,'XTickLabel',ctext)
% print('-dpng',[figp 'corr_rmse_bias_clims_200_obsglm_neg11_4rt_together.png'])

%% Tables
cor = squeeze(skill(1,:,:));
Tcorr = array2table(cor,'RowNames',simtext,'VariableNames',ctext);

rmse = squeeze(skill(2,:,:));
Trmse = array2table(rmse,'RowNames',simtext,'VariableNames',ctext);

nstd = squeeze(skill(8,:,:));
Tnstd = array2table(nstd,'RowNames',simtext,'VariableNames',ctext);

urmse = squeeze(skill(9,:,:));
Turmse = array2table(urmse,'RowNames',simtext,'VariableNames',ctext);

trmse = squeeze(skill(10,:,:));
Ttrmse = array2table(trmse,'RowNames',simtext,'VariableNames',ctext);

bias = squeeze(skill(11,:,:));
Tbias = array2table(bias,'RowNames',simtext,'VariableNames',ctext);

writetable(Tcorr,[sfile 'corr_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);
writetable(Trmse,[sfile 'rmse_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);
writetable(Tnstd,[sfile 'nstd_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);
writetable(Turmse,[sfile 'urmse_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);
writetable(Ttrmse,[sfile 'trmse_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);
writetable(Tbias,[sfile 'bias_hist_clims_200_obsglm100_neg11_4rt.csv'],'WriteRowNames',true);

