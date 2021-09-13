% Linear regressions of zmeso biomass with surf chl 
% Tropical regions, annual climatology
% Historic and SSP585

clear all
close all

%%
load('Hist_coeffs_mod_chl_Ntemp_SprSum_200.mat');
Hist.mod = mod;
Hist.slop = slop;
Hist.int = int;

clear mod slop int

%%
load('SSP585_coeffs_mod_chl_Ntemp_SprSum_200.mat');
SSP.mod = mod;
SSP.slop = slop;
SSP.int = int;

clear mod slop int

%% predict mesozoo over chl values
% Models have a very large range
% Don't know units of obs
chl = logspace(-5,-2);
hzoo = Hist.int + Hist.slop * log10(chl);
szoo = SSP.int + SSP.slop * log10(chl);

%% color same as Taylor diagrams
cm=[0 0.7 0;...   %g
    0 0 0.65;...  %b
    0.7 0 1;...   %purple
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.35 0.35 0.35]; %grey
set(groot,'defaultAxesColorOrder',cm);

%% figures
figure(1)
subplot(2,2,1)
plot(log10(chl),hzoo(1:5,:),'LineWidth',2)
legend(SSP.mod)
legend('location','northwest')
xlim([-5 -2])
xlabel('log_1_0 chl (kg m^-^3)')
ylabel('log_1_0 zmeso (mgC m^-^2)')
title('Hist')

subplot(2,2,3)
plot(log10(chl),hzoo(1:5,:),'LineWidth',2)
% legend(SSP.mod)
% legend('location','northwest')
xlim([-5 -2])
ylim([-5 10])
xlabel('log_1_0 chl (kg m^-^3)')
ylabel('log_1_0 zmeso (mgC m^-^2)')
title('Hist')

subplot(2,2,2)
plot(log10(chl),szoo,'LineWidth',2)
xlim([-5 -2])
xlabel('log_1_0 chl (kg m^-^3)')
ylabel('log_1_0 zmeso (mgC m^-^2)')
title('SSP 585')

subplot(2,2,4)
plot(log10(chl),szoo,'LineWidth',2)
xlim([-5 -2])
ylim([-5 10])
xlabel('log_1_0 chl (kg m^-^3)')
ylabel('log_1_0 zmeso (mgC m^-^2)')
title('SSP 585')
print('-dpng','Hist_SSP585_lm_zm_chl_Ntemp_SprSum_clim_200.png')

%% linear regression of hist vs. future
% MISSING TOOLBOX?
mdl = fitlm(Hist.int(1:5),SSP.int);
ypred = predict(mdl,Hist.int(1:5));

mdl2 = fitlm(Hist.int(2:5),SSP.int(2:5));
ypred2 = predict(mdl2,Hist.int(2:5));

test=mdl.Coefficients;
tt=table2array(test);

%%
figure(2)
subplot(2,2,1)
for i=1:5
    plot(Hist.int(i),SSP.int(i),'.','MarkerSize',20); hold on;
end
%plot(Hist.int(1:5),ypred,'--k'); hold on;
legend(SSP.mod)
legend('location','northwest')
%xlim([-5 -2])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('Zmeso vs. chl')
text(40,10,'R^2 = 0.99');
text(40,5,'p = 5.18e-5');

subplot(2,2,3)
for i=1:5
    plot(Hist.int(i),SSP.int(i),'.','MarkerSize',20); hold on;
end
%plot(Hist.int(2:5),ypred2,'--k'); hold on;
% legend(SSP.mod)
% legend('location','southeast')
xlim([3 14])
xlabel('Historic chl sensitivity')
ylabel('Future chl sensitivity')
title('Zmeso vs. chl')
text(10,4,'R^2 = 0.86');
text(10,3,'p = 0.07');
print('-dpng','Hist_SSP585_sens_zm_chl_Ntemp_SprSum_clim_200.png')


