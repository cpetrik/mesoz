% CMIP6 model Q10 values

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cqpg = 1.719;
cqzg = 1.0;
cqpm = 1.0;
cqzm = 1.0;
cqdr = 2.168;

%% CMCC
mqpg = 2;
mqzg = 2;
mqpm = 2;
mqzm = 2;
mqdr = 2;

%% PISCES
pqpg = 1.895;
pqzg = 2.139;
pqpm = 1.0;
pqzm = 2.139;
pqdr = 1.9;

%% GFDL
gqpg = 1.878;
gqzg = 1.878;
gqpm = 1.0;
gqzm = 1.878;
gqdr = 1.878;

%% UKESM
uqpg = 1.895;
uqzg = 1.0;
uqpm = 1.0;
uqzm = 1.0;
uqdr = 1.895;

%% Plot ratio of rates as fn of diff in temp
temp1 = 0:0.75:10;
temp2 = -0.2:0.75:10;
temp3 = -0.4:0.75:10;
temp4 = -0.6:0.75:10;
temp5 = -0.8:0.75:10;

% phyto growth
cpg = cqpg .^ (temp5./10);
mpg = mqpg .^ (temp2./10);
ppg = pqpg .^ (temp3./10);
gpg = gqpg .^ (temp4./10);
upg = uqpg .^ (temp1./10);

% zoo grazing
czg = cqzg .^ (temp2./10);
mzg = mqzg .^ (temp3./10);
pzg = pqzg .^ (temp4./10);
gzg = gqzg .^ (temp5./10);
uzg = uqzg .^ (temp1./10);

% phyto mort
cpm = cqpm .^ (temp3./10);
mpm = mqpm .^ (temp4./10);
ppm = pqpm .^ (temp5./10);
gpm = gqpm .^ (temp1./10);
upm = uqpm .^ (temp2./10);

% zoo mort
czm = cqzm .^ (temp4./10);
mzm = mqzm .^ (temp5./10);
pzm = pqzm .^ (temp1./10);
gzm = gqzm .^ (temp2./10);
uzm = uqzm .^ (temp3./10);

% detritus remin
cdr = cqdr .^ (temp5./10);
mdr = mqdr .^ (temp1./10);
pdr = pqdr .^ (temp2./10);
gdr = gqdr .^ (temp3./10);
udr = uqdr .^ (temp4./10);

%% Plots 
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...   %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0];                     %black
cb=[34/255 136/255 51/255;...   %green
    153/255 153/255 51/255;...   %olive
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255;...    %blue
    170/255 51/255 119/255];    %purple

set(groot,'defaultAxesColorOrder',cb);

%%
%figure(1)
figure('Units','inches','Position',[1 2 6 9]);
subplot(2,1,1)
plot(temp2,czg,'--*','MarkerSize',5); hold on;
plot(temp3,mzg,'--*','MarkerSize',5); hold on;
plot(temp4,pzg,'--*','MarkerSize',5); hold on;
plot(temp5,gzg,'--*','MarkerSize',5); hold on;
plot(temp1,uzg,'--*','MarkerSize',5); hold on;

plot(temp5,cpg,'--.','MarkerSize',12); hold on;
plot(temp2,mpg,'--.','MarkerSize',12); hold on;
plot(temp3,ppg,'--.','MarkerSize',12); hold on;
plot(temp4,gpg,'--.','MarkerSize',12); hold on;
plot(temp1,upg,'--.','MarkerSize',12); hold on;
legend({'zCAN','zCMCC','zCNRM/IPSL','zGFDL','zUK',...
    'pCAN','pCMCC','pCNRM/IPSL','pGFDL','pUK'},'Location','northwest','NumColumns',2);
ylim([0.8 2.2])
xlim([0 10])
xlabel('Difference in temperature')
ylabel('Ratio of rates')
title('Growth or Grazing')

subplot(2,1,2)
plot(temp5,cdr,'--v','MarkerSize',5); hold on;
plot(temp1,mdr,'--v','MarkerSize',5); hold on;
plot(temp2,pdr,'--v','MarkerSize',5); hold on;
plot(temp3,gdr,'--v','MarkerSize',5); hold on;
plot(temp4,udr,'--v','MarkerSize',5); hold on;
lgd = legend({'dCAN','dCMCC','dCNRM/IPSL','dGFDL','dUK'},...
    'location','northwest');
lgd.AutoUpdate = 'off';

plot(temp4,czm,'--*','MarkerSize',5); hold on;
plot(temp5,mzm,'--*','MarkerSize',5); hold on;
plot(temp1,pzm,'--*','MarkerSize',5); hold on;
plot(temp2,gzm,'--*','MarkerSize',5); hold on;
plot(temp3,uzm,'--*','MarkerSize',5); hold on;

plot(temp3,cpm,'--.','MarkerSize',12); hold on;
plot(temp4,mpm,'--.','MarkerSize',12); hold on;
plot(temp5,ppm,'--.','MarkerSize',12); hold on;
plot(temp1,gpm,'--.','MarkerSize',12); hold on;
plot(temp2,upm,'--.','MarkerSize',12); hold on;
ylim([0.8 2.2])
xlim([0 10])
xlabel('Difference in temperature')
ylabel('Ratio of rates')
title('Mortality or Remineralization')
print('-dpng',[ppath 'temp_dep_q10_bgc.png'])

%%
%figure(1)
figure('Units','inches','Position',[1 2 13 6]);
subplot(1,2,1)
plot(temp2,czg,'--*','MarkerSize',5); hold on;
plot(temp3,mzg,'--*','MarkerSize',5); hold on;
plot(temp4,pzg,'--*','MarkerSize',5); hold on;
plot(temp5,gzg,'--*','MarkerSize',5); hold on;
plot(temp1,uzg,'--*','MarkerSize',5); hold on;

plot(temp5,cpg,'--.','MarkerSize',12); hold on;
plot(temp2,mpg,'--.','MarkerSize',12); hold on;
plot(temp3,ppg,'--.','MarkerSize',12); hold on;
plot(temp4,gpg,'--.','MarkerSize',12); hold on;
plot(temp1,upg,'--.','MarkerSize',12); hold on;
legend({'zCAN','zCMCC','zCNRM/IPSL','zGFDL','zUK',...
    'pCAN','pCMCC','pCNRM/IPSL','pGFDL','pUK'},'Location','northwest','NumColumns',2);
ylim([0.8 2.2])
xlim([0 10])
xlabel('Difference in temperature')
ylabel('Ratio of rates')
title('Growth or Grazing')

subplot(1,2,2)
plot(temp5,cdr,'--v','MarkerSize',5); hold on;
plot(temp1,mdr,'--v','MarkerSize',5); hold on;
plot(temp2,pdr,'--v','MarkerSize',5); hold on;
plot(temp3,gdr,'--v','MarkerSize',5); hold on;
plot(temp4,udr,'--v','MarkerSize',5); hold on;
lgd = legend({'dCAN','dCMCC','dCNRM/IPSL','dGFDL','dUK'},...
    'location','northwest');
lgd.AutoUpdate = 'off';

plot(temp4,czm,'--*','MarkerSize',5); hold on;
plot(temp5,mzm,'--*','MarkerSize',5); hold on;
plot(temp1,pzm,'--*','MarkerSize',5); hold on;
plot(temp2,gzm,'--*','MarkerSize',5); hold on;
plot(temp3,uzm,'--*','MarkerSize',5); hold on;

plot(temp3,cpm,'--.','MarkerSize',12); hold on;
plot(temp4,mpm,'--.','MarkerSize',12); hold on;
plot(temp5,ppm,'--.','MarkerSize',12); hold on;
plot(temp1,gpm,'--.','MarkerSize',12); hold on;
plot(temp2,upm,'--.','MarkerSize',12); hold on;
ylim([0.8 2.2])
xlim([0 10])
xlabel('Difference in temperature')
ylabel('Ratio of rates')
title('Mortality or Remineralization')
print('-dpng',[ppath 'temp_dep_q10_bgc_v2.png'])


