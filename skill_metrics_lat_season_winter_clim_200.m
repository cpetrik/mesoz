% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

%% 
load('obs_mod_chl_spol_winter_clim_200.mat')
obs = obsmodchl(:,3);
model = obsmodchl(:,4:8);

%% Skill metric = weighted sum of squares for multivariate

metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

skill=NaN*ones(14,5);
for j=1:5
    o=obs;
    p=model(:,j);
    
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
    skill(1,j) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    skill(2,j) = sqrt(num/n);
    
    % reliability index  %%What about predicted zero values?
    %q1=nansum((log(o./p)).^2);
    q1=nansum((o-p).^2);
    skill(3,j) = exp(sqrt(q1/n));
    
    % average error
    skill(4,j) = nansum(p-o) / n;
    
    % average absolute error
    skill(5,j) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    skill(6,j) = (num1-num2)/num1;
    
    % Taylor R
    num=nansum((o-omean).*(p-pmean));
    skill(7,j) = num/(n*osig*psig);
    
    % Taylor normalized std
    skill(8,j) = psig/osig;
    
    % unbiased root mean square difference
    % sign tells difference between model and obs std
    q1=nansum(((p-pmean)-(o-omean)).^2);
    if (psig>=osig)
        skill(9,j) = sqrt(q1/n);
    else
        skill(9,j) = -1*sqrt(q1/n);
    end
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    skill(10,j) = sqrt(q1/n);
    
    % normalized bias
    skill(11,j) = (pmean(1,1)-omean(1,1));%./osig;
    
    % Joliff et al 2009 S1
    R=skill(7,j);
    s=skill(8,j);
    num=2*(1+R);
    den=(s + (1/s)).^2;
    skill(12,j) = 1 - (num/den);
    
    % Joliff et al 2009 S2
    num=(1+R).^4;
    den=4*(s + (1/s)).^2;
    skill(13,j) = 1 - (num/den);
    
    % Joliff et al 2009 S3
    q1=exp(-((s-1).^2) / 0.18);
    q2=(1+R)/2;
    skill(14,j) = 1 - (q1*q2);
    
end

%% Results
simtex = {'CAN','CNRM','GFDL','IPSL','UK'};

% Bar graphs
figure(1)
subplot(3,2,1)
bar(skill(1,:),'k')
set(gca,'XTickLabel',simtex)

subplot(3,2,3)
bar(skill(2,:),'k')
set(gca,'XTickLabel',simtex)

subplot(3,2,5)
bar(skill(6,:),'k')
set(gca,'XTickLabel',simtex)
ylim([-35 5])

subplot(3,2,2)
bar(skill(1,:),'k')
text(-1,0.15,'Winter Correlation coefficient','FontWeight','Bold','HorizontalAlignment','center')
set(gca,'XTickLabel',simtex)
ylim([-0.15 0.1])

subplot(3,2,4)
bar(skill(2,:),'k')
text(-1,1.0,'Root mean square error','FontWeight','Bold','HorizontalAlignment','center')
set(gca,'XTickLabel',simtex)
ylim([0 0.9])

subplot(3,2,6)
bar(skill(6,:),'k')
text(-1,0.4,'Modeling efficiency','FontWeight','Bold','HorizontalAlignment','center')
set(gca,'XTickLabel',simtex)
ylim([-0.75 0.25])
stamp('S Polar Winter')
print('-dpng','corr_rmse_mef_spol_winter_clim_200.png')


%% Taylor diagram 
% cm={[1 0.5 0],...   %orange
%     [0.5 0.5 0],... %tan/army
%     [0 0.7 0],...   %g
%     [0 1 1],...     %c
%     [0 0 0.75],...  %b
%     [0.5 0 1],...   %purple
%     [1 0 1],...     %m
%     [1 0 0],...     %r
%     [0.5 0 0],...   %maroon
%     [0.35 0.35 0.35]}; %grey


[rmsd,it]=sort(skill(10,:),'descend');  %tot RMSD
theta=acos(skill(7,:));                 %Taylor R
rho=skill(8,:);                         %Taylor norm std
simtex{6}='obs';

cm={[0 0.7 0],...   %g
    [0 0 0.65],...  %b
    [0.7 0 1],...   %purple
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.35 0.35 0.35]}; %grey

tr=0;
rr=1;

%% Taylor diagram using corr coeff calculated
theta2=acos(skill(1,:));

figure(3)
subplot(1,2,1)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 5])
title('Annual climatology')

subplot(1,2,2)
h0=polarplot(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:5
    h=polarplot(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polarplot(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 90 0 1.25])
title('Annual climatology')
legend([' ' simtex])
legend('location','northeast')
print('-dpng','Taylor_spol_winter_clim_200_clim_200')


%% 









