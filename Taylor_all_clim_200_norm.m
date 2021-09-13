% Calculate different skill metrics for each ESM
% log transform biomass

clear all
close all

%%
load('skill_model_obsglm_clim_zmeso200.mat')

%% Standardization
% 1st take log
Acomb = log(comb(:,3:8) + eps);
Dcomb = log(dcomb(:,3:8) + eps);
Jcomb = log(jcomb(:,3:8) + eps);
Mcomb = log(mcomb(:,3:8) + eps);
Scomb = log(scomb(:,3:8) + eps);
% then stdize by (x-mean)/std
%Acomb = normalize(Acomb); - didn't work
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

%% Calculate statistics for Taylor diagram
% The first array element corresponds to the reference series for the
% while the second is that for the predicted series.

%Annual means
for m=1:5
    ref=Acomb(:,6);
    pred=Acomb(:,m);
    pred(isnan(ref))=[];
    ref(isnan(ref))=[];
    ref(isnan(pred))=[];
    pred(isnan(pred))=[];
    eval(['atay_st' num2str(m) ' = taylor_statistics(pred,ref);']);
end

%MAM means
for m=1:5
    ref=Mcomb(:,6);
    pred=Mcomb(:,m);
    pred(isnan(ref))=[];
    ref(isnan(ref))=[];
    ref(isnan(pred))=[];
    pred(isnan(pred))=[];
    eval(['mtay_st' num2str(m) ' = taylor_statistics(pred,ref);']);
end

%JJA means
for m=1:5
    ref=Jcomb(:,6);
    pred=Jcomb(:,m);
    pred(isnan(ref))=[];
    ref(isnan(ref))=[];
    ref(isnan(pred))=[];
    pred(isnan(pred))=[];
    eval(['jtay_st' num2str(m) ' = taylor_statistics(pred,ref);']);
end

%SON means
for m=1:5
    ref=Scomb(:,6);
    pred=Scomb(:,m);
    pred(isnan(ref))=[];
    ref(isnan(ref))=[];
    ref(isnan(pred))=[];
    pred(isnan(pred))=[];
    eval(['stay_st' num2str(m) ' = taylor_statistics(pred,ref);']);
end

%DJF means
for m=1:5
    ref=Dcomb(:,6);
    pred=Dcomb(:,m);
    pred(isnan(ref))=[];
    ref(isnan(ref))=[];
    ref(isnan(pred))=[];
    pred(isnan(pred))=[];
    eval(['dtay_st' num2str(m) ' = taylor_statistics(pred,ref);']);
end

%% Store statistics in arrays
Asdev = [atay_st1.sdev(1); atay_st1.sdev(2); ...
    atay_st2.sdev(2); atay_st3.sdev(2); atay_st4.sdev(2); atay_st5.sdev(2)];
Acrmsd = [atay_st1.crmsd(1); atay_st1.crmsd(2); ...
    atay_st2.crmsd(2); atay_st3.crmsd(2); atay_st4.crmsd(2); atay_st5.crmsd(2)];
Accoef = [atay_st1.ccoef(1); atay_st1.ccoef(2); ...
    atay_st2.ccoef(2); atay_st3.ccoef(2); atay_st4.ccoef(2); atay_st5.ccoef(2)];

Msdev = [mtay_st1.sdev(1); mtay_st1.sdev(2); ...
    mtay_st2.sdev(2); mtay_st3.sdev(2); mtay_st4.sdev(2); mtay_st5.sdev(2)];
Mcrmsd = [mtay_st1.crmsd(1); mtay_st1.crmsd(2); ...
    mtay_st2.crmsd(2); mtay_st3.crmsd(2); mtay_st4.crmsd(2); mtay_st5.crmsd(2)];
Mccoef = [mtay_st1.ccoef(1); mtay_st1.ccoef(2); ...
    mtay_st2.ccoef(2); mtay_st3.ccoef(2); mtay_st4.ccoef(2); mtay_st5.ccoef(2)];

Jsdev = [jtay_st1.sdev(1); jtay_st1.sdev(2); ...
    jtay_st2.sdev(2); jtay_st3.sdev(2); jtay_st4.sdev(2); jtay_st5.sdev(2)];
Jcrmsd = [jtay_st1.crmsd(1); jtay_st1.crmsd(2); ...
    jtay_st2.crmsd(2); jtay_st3.crmsd(2); jtay_st4.crmsd(2); jtay_st5.crmsd(2)];
Jccoef = [jtay_st1.ccoef(1); jtay_st1.ccoef(2); ...
    jtay_st2.ccoef(2); jtay_st3.ccoef(2); jtay_st4.ccoef(2); jtay_st5.ccoef(2)];

Ssdev = [stay_st1.sdev(1); stay_st1.sdev(2); ...
    stay_st2.sdev(2); stay_st3.sdev(2); stay_st4.sdev(2); stay_st5.sdev(2)];
Scrmsd = [stay_st1.crmsd(1); stay_st1.crmsd(2); ...
    stay_st2.crmsd(2); stay_st3.crmsd(2); stay_st4.crmsd(2); stay_st5.crmsd(2)];
Sccoef = [stay_st1.ccoef(1); stay_st1.ccoef(2); ...
    stay_st2.ccoef(2); stay_st3.ccoef(2); stay_st4.ccoef(2); stay_st5.ccoef(2)];

Dsdev = [dtay_st1.sdev(1); dtay_st1.sdev(2); ...
    dtay_st2.sdev(2); dtay_st3.sdev(2); dtay_st4.sdev(2); dtay_st5.sdev(2)];
Dcrmsd = [dtay_st1.crmsd(1); dtay_st1.crmsd(2); ...
    dtay_st2.crmsd(2); dtay_st3.crmsd(2); dtay_st4.crmsd(2); dtay_st5.crmsd(2)];
Dccoef = [dtay_st1.ccoef(1); dtay_st1.ccoef(2); ...
    dtay_st2.ccoef(2); dtay_st3.ccoef(2); dtay_st4.ccoef(2); dtay_st5.ccoef(2)];

%% Produce the Taylor diagram.
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35];
%
label = ['ref' simtex];
[hp, ht, axl] = taylor_diagram(Asdev,Acrmsd,Accoef,'numberPanels', 2, ...
    'markerLabel',label, 'colormap', cm, ...
    'tickRMS',0.0:10.0:80.0,'tickRMSangle',150.0, ...
    'colRMS','m', 'styleRMS', ':', 'widthRMS', 2.0, 'titleRMS', 'off', ...
    'tickSTD',0.0:20.0:60.0, 'limSTD',60.0, ...
    'colSTD','b', 'styleSTD', '-.', 'widthSTD', 1.0, ...
    'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0);
% Write plot to file
writepng(gcf,'taylor_annual_norm.png');

%% Results
figp ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';
cm=[0 0.7 0;...   %g
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey
set(groot,'defaultAxesColorOrder',cm);



