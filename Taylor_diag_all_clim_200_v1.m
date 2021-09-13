% Calculate different skill metrics for each ESM
% log transform biomass
% Taylor diagram code

clear all
close all

%%
load('skill_all_clim_200.mat');

%%
id1 = find(~isnan(model(:,1)));
id2 = find(~isnan(model(:,2)));
id3 = find(~isnan(model(:,3)));
id4 = find(~isnan(model(:,4)));
id5 = find(~isnan(model(:,5)));

%% Calculate statistics for Taylor diagram
% The first array element corresponds to the reference series for the
% while the second is that for the predicted series.
taylor_stats1 = taylor_statistics(model(id1,1),obs(id1));
taylor_stats2 = taylor_statistics(model(id2,2),obs(id2));
taylor_stats3 = taylor_statistics(model(id3,3),obs(id3));
taylor_stats4 = taylor_statistics(model(id4,4),obs(id4));
taylor_stats5 = taylor_statistics(model(id5,5),obs(id5));

%% Store statistics in arrays
sdev = [taylor_stats1.sdev(1); taylor_stats1.sdev(2); ...
    taylor_stats2.sdev(2); taylor_stats3.sdev(2);...
    taylor_stats4.sdev(2); taylor_stats5.sdev(2)];
crmsd = [taylor_stats1.crmsd(1); taylor_stats1.crmsd(2); ...
    taylor_stats2.crmsd(2); taylor_stats3.crmsd(2);...
    taylor_stats4.crmsd(2); taylor_stats5.crmsd(2)];
ccoef = [taylor_stats1.ccoef(1); taylor_stats1.ccoef(2); ...
    taylor_stats2.ccoef(2); taylor_stats3.ccoef(2);...
    taylor_stats4.ccoef(2); taylor_stats5.ccoef(2)];

%%
label = {'Observation', 'CAN', 'CNRM', 'GFDL', 'IPSL', 'UK'};

%% Produce the Taylor diagram.
%
% Note that the first index corresponds to the reference series for the
% diagram. For example sdev(1) is the standard deviation of the reference
% series and sdev(2:4) are the standard deviations of the other series.
% The value of sdev(1) is used to define the origin of the RMSD contours.
% The other values are used to plot the points (total of 3) that appear in
% the diagram.

%[hp, ht, axl] = taylor_diagram(sdev,crmsd,ccoef);

[hp, ht, axl] = taylor_diagram(sdev,crmsd,ccoef, ...
    'markerLabel',label, 'markerColor', 'r', 'markerLegend', 'on', ...
    'colRMS','m', 'styleRMS', ':', 'widthRMS', 2.0, 'titleRMS', 'on', ...
    'colSTD','b', 'styleSTD', '-.', 'widthSTD', 1.0, 'titleSTD', 'on', ...
    'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0, 'titleCOR', 'on', ...
    'markerSize',14, 'alpha', 0.0);

% [hp, ht, axl] = taylor_diagram(sdev,crmsd,ccoef, ...
%     'markerLabel',label, 'markerColor', 'r', 'markerLegend', 'on', ...
%     'tickRMS',0.0:5.0, ...
%     'colRMS','m', 'styleRMS', ':', 'widthRMS', 2.0, 'titleRMS', 'on', ...
%     'tickSTD',0.0:5.0, 'limSTD',5.0, ...
%     'colSTD','b', 'styleSTD', '-.', 'widthSTD', 1.0, 'titleSTD', 'on', ...
%     'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0, 'titleCOR', 'on', ...
%     'markerSize',14, 'alpha', 0.0);

% Write plot to file
%writepng(gcf,'taylor1.png');

%% Use other function

% rmsd=(skill(10,:));       %tot RMSD
% theta=acos(skill(7,:));   %Taylor R
% rho=skill(8,:);           %Stdev

%[hp ht axl] = taylordiag(STDs,RMSs,CORs,['option',value])
[hp ht axl] = taylordiag(skill(8,:),skill(10,:),skill(7,:))



