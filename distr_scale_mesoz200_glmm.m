% CMIP6 output 
% mesoz 200m integrations
% Histograms of ESM with diff transformations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_space_means_50yr_zmeso200_same_orientation.mat');

%% Log-trans add eps
cz_e = log(cz+eps);
mz_e = log(mz+eps);
nz_e = log(nz+eps);
gz_e = log(gz+eps);
iz_e = log(iz+eps);
uz_e = log(uz+eps);
oz_e = log(oz+eps);

%% Log-trans add 1/2 min
cmin = 0.5 * min(cz(cz(:) > 0));
mmin = 0.5 * min(mz(mz(:) > 0));
nmin = 0.5 * min(nz(nz(:) > 0));
gmin = 0.5 * min(gz(gz(:) > 0));
imin = 0.5 * min(iz(iz(:) > 0));
umin = 0.5 * min(uz(uz(:) > 0));
omin = 0.5 * min(oz(oz(:) > 0));

cz_h = log(cz+cmin);
mz_h = log(mz+mmin);
nz_h = log(nz+nmin);
gz_h = log(gz+gmin);
iz_h = log(iz+imin);
uz_h = log(uz+umin);
oz_h = log(oz+omin);

%% Log-trans add 1E-8
cz_8 = log(cz+1e-8);
mz_8 = log(mz+1e-8);
nz_8 = log(nz+1e-8);
gz_8 = log(gz+1e-8);
iz_8 = log(iz+1e-8);
uz_8 = log(uz+1e-8);
oz_8 = log(oz+1e-8);

%% Log-trans add 1
cz_1 = log(cz+1);
mz_1 = log(mz+1);
nz_1 = log(nz+1);
gz_1 = log(gz+1);
iz_1 = log(iz+1);
uz_1 = log(uz+1);
oz_1 = log(oz+1);

%% sqrt-trans
cz_sq = sqrt(cz);
mz_sq = sqrt(mz);
nz_sq = sqrt(nz);
gz_sq = sqrt(gz);
iz_sq = sqrt(iz);
uz_sq = sqrt(uz);
oz_sq = sqrt(oz);

%% fourth-root-trans
cz_4th = (cz).^(1/4);
mz_4th = (mz).^(1/4);
nz_4th = (nz).^(1/4);
gz_4th = (gz).^(1/4);
iz_4th = (iz).^(1/4);
uz_4th = (uz).^(1/4);
oz_4th = (oz).^(1/4);

%% Histograms
% Untrans
figure(1)
subplot(3,3,1)
histogram(cz(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz(:),'FaceColor','k');
title({'Untrans' , 'CMCC'})

subplot(3,3,3)
histogram(nz(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_untrans.png'])

%% Log + eps
figure(2)
subplot(3,3,1)
histogram(cz_e(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_e(:),'FaceColor','k');
title({'log(b+eps)' , 'CMCC'})

subplot(3,3,3)
histogram(nz_e(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_e(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_e(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_e(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_e(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_logeps.png'])

%% Log + 1/2 min
figure(3)
subplot(3,3,1)
histogram(cz_h(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_h(:),'FaceColor','k');
title({'log(b+0.5(min))' , 'CMCC'})

subplot(3,3,3)
histogram(nz_h(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_h(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_h(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_h(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_h(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_loghalfmin.png'])

%% Log + 1e-8
figure(4)
subplot(3,3,1)
histogram(cz_8(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_8(:),'FaceColor','k');
title({'log(b+1e-8)' , 'CMCC'})

subplot(3,3,3)
histogram(nz_8(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_8(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_8(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_8(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_8(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_log1e-8.png'])

%% Log + 1
figure(5)
subplot(3,3,1)
histogram(cz_1(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_1(:),'FaceColor','k');
title({'log(b+1)' , 'CMCC'})

subplot(3,3,3)
histogram(nz_1(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_1(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_1(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_1(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_1(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_log1.png'])

%% Sqrt
figure(6)
subplot(3,3,1)
histogram(cz_sq(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_sq(:),'FaceColor','k');
title({'Sqrt' , 'CMCC'})

subplot(3,3,3)
histogram(nz_sq(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_sq(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_sq(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_sq(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_sq(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_sqrt.png'])

%% 4th-root
figure(7)
subplot(3,3,1)
histogram(cz_4th(:),'FaceColor','k');
title('CAN')

subplot(3,3,2)
histogram(mz_4th(:),'FaceColor','k');
title({'4th root' , 'CMCC'})

subplot(3,3,3)
histogram(nz_4th(:),'FaceColor','k');
title('CNRM')

subplot(3,3,4)
histogram(gz_4th(:),'FaceColor','k');
title('GFDL')

subplot(3,3,5)
histogram(iz_4th(:),'FaceColor','k');
title('IPSL')

subplot(3,3,6)
histogram(uz_4th(:),'FaceColor','k');
title('UK')

subplot(3,3,8)
histogram(oz_4th(:),'FaceColor','k');
title('obsGLMM')
print('-dpng',[ppath 'Distr_all_clim_obsglm_4rt.png'])

%% tests
[h(1,1),p(1,1),k(1,1),c(1,1)] = kstest(cz(:));
[h(1,2),p(1,2),k(1,2),c(1,2)] = kstest(cz_e(:));
[h(1,3),p(1,3),k(1,3),c(1,3)] = kstest(cz_h(:));
[h(1,4),p(1,4),k(1,4),c(1,4)] = kstest(cz_8(:));
[h(1,5),p(1,5),k(1,5),c(1,5)] = kstest(cz_1(:));
[h(1,6),p(1,6),k(1,6),c(1,6)] = kstest(cz_sq(:));
[h(1,7),p(1,7),k(1,7),c(1,7)] = kstest(cz_4th(:));

[h(2,1),p(2,1),k(2,1),c(2,1)] = kstest(mz(:));
[h(2,2),p(2,2),k(2,2),c(2,2)] = kstest(mz_e(:));
[h(2,3),p(2,3),k(2,3),c(2,3)] = kstest(mz_h(:));
[h(2,4),p(2,4),k(2,4),c(2,4)] = kstest(mz_8(:));
[h(2,5),p(2,5),k(2,5),c(2,5)] = kstest(mz_1(:));
[h(2,6),p(2,6),k(2,6),c(2,6)] = kstest(mz_sq(:));
[h(2,7),p(2,7),k(2,7),c(2,7)] = kstest(mz_4th(:));

[h(3,1),p(3,1),k(3,1),c(3,1)] = kstest(nz(:));
[h(3,2),p(3,2),k(3,2),c(3,2)] = kstest(nz_e(:));
[h(3,3),p(3,3),k(3,3),c(3,3)] = kstest(nz_h(:));
[h(3,4),p(3,4),k(3,4),c(3,4)] = kstest(nz_8(:));
[h(3,5),p(3,5),k(3,5),c(3,5)] = kstest(nz_1(:));
[h(3,6),p(3,6),k(3,6),c(3,6)] = kstest(nz_sq(:));
[h(3,7),p(3,7),k(3,7),c(3,7)] = kstest(nz_4th(:));

[h(4,1),p(4,1),k(4,1),c(4,1)] = kstest(gz(:));
[h(4,2),p(4,2),k(4,2),c(4,2)] = kstest(gz_e(:));
[h(4,3),p(4,3),k(4,3),c(4,3)] = kstest(gz_h(:));
[h(4,4),p(4,4),k(4,4),c(4,4)] = kstest(gz_8(:));
[h(4,5),p(4,5),k(4,5),c(4,5)] = kstest(gz_1(:));
[h(4,6),p(4,6),k(4,6),c(4,6)] = kstest(gz_sq(:));
[h(4,7),p(4,7),k(4,7),c(4,7)] = kstest(gz_4th(:));

[h(5,1),p(5,1),k(5,1),c(5,1)] = kstest(iz(:));
[h(5,2),p(5,2),k(5,2),c(5,2)] = kstest(iz_e(:));
[h(5,3),p(5,3),k(5,3),c(5,3)] = kstest(iz_h(:));
[h(5,4),p(5,4),k(5,4),c(5,4)] = kstest(iz_8(:));
[h(5,5),p(5,5),k(5,5),c(5,5)] = kstest(iz_1(:));
[h(5,6),p(5,6),k(5,6),c(5,6)] = kstest(iz_sq(:));
[h(5,7),p(5,7),k(5,7),c(5,7)] = kstest(iz_4th(:));

[h(6,1),p(6,1),k(6,1),c(6,1)] = kstest(uz(:));
[h(6,2),p(6,2),k(6,2),c(6,2)] = kstest(uz_e(:));
[h(6,3),p(6,3),k(6,3),c(6,3)] = kstest(uz_h(:));
[h(6,4),p(6,4),k(6,4),c(6,4)] = kstest(uz_8(:));
[h(6,5),p(6,5),k(6,5),c(6,5)] = kstest(uz_1(:));
[h(6,6),p(6,6),k(6,6),c(6,6)] = kstest(uz_sq(:));
[h(6,7),p(6,7),k(6,7),c(6,7)] = kstest(uz_4th(:));

[h(7,1),p(7,1),k(7,1),c(7,1)] = kstest(oz(:));
[h(7,2),p(7,2),k(7,2),c(7,2)] = kstest(oz_e(:));
[h(7,3),p(7,3),k(7,3),c(7,3)] = kstest(oz_h(:));
[h(7,4),p(7,4),k(7,4),c(7,4)] = kstest(oz_8(:));
[h(7,5),p(7,5),k(7,5),c(7,5)] = kstest(oz_1(:));
[h(7,6),p(7,6),k(7,6),c(7,6)] = kstest(oz_sq(:));
[h(7,7),p(7,7),k(7,7),c(7,7)] = kstest(oz_4th(:));

%% Min diff betwen k and c
dkc = k-c;
mdkc = mean(dkc);
ndkc = median(dkc);

%no CAN
mdkc2 = mean(dkc(2:7,:));
ndkc2 = median(dkc(2:7,:));

%% Normplots
close all
% Untrans
figure(1)
subplot(3,3,1)
normplot(cz(:));
title('CAN')

subplot(3,3,2)
normplot(mz(:));
title({'Untrans' , 'CMCC'})

subplot(3,3,3)
normplot(nz(:));
title('CNRM')

subplot(3,3,4)
normplot(gz(:));
title('GFDL')

subplot(3,3,5)
normplot(iz(:));
title('IPSL')

subplot(3,3,6)
normplot(uz(:));
title('UK')

subplot(3,3,8)
normplot(oz(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_untrans.png'])

%% Log + eps
figure(2)
subplot(3,3,1)
normplot(cz_e(:));
title('CAN')

subplot(3,3,2)
normplot(mz_e(:));
title({'log(b+eps)' , 'CMCC'})

subplot(3,3,3)
normplot(nz_e(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_e(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_e(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_e(:));
title('UK')

subplot(3,3,8)
normplot(oz_e(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_logeps.png'])

%% Log + 1/2 min
figure(3)
subplot(3,3,1)
normplot(cz_h(:));
title('CAN')

subplot(3,3,2)
normplot(mz_h(:));
title({'log(b+0.5(min))' , 'CMCC'})

subplot(3,3,3)
normplot(nz_h(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_h(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_h(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_h(:));
title('UK')

subplot(3,3,8)
normplot(oz_h(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_loghalfmin.png'])

%% Log + 1e-8
figure(4)
subplot(3,3,1)
normplot(cz_8(:));
title('CAN')

subplot(3,3,2)
normplot(mz_8(:));
title({'log(b+1e-8)' , 'CMCC'})

subplot(3,3,3)
normplot(nz_8(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_8(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_8(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_8(:));
title('UK')

subplot(3,3,8)
normplot(oz_8(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_log1e-8.png'])

%% Log + 1
figure(5)
subplot(3,3,1)
normplot(cz_1(:));
title('CAN')

subplot(3,3,2)
normplot(mz_1(:));
title({'log(b+1)' , 'CMCC'})

subplot(3,3,3)
normplot(nz_1(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_1(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_1(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_1(:));
title('UK')

subplot(3,3,8)
normplot(oz_1(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_log1.png'])

%% Sqrt
figure(6)
subplot(3,3,1)
normplot(cz_sq(:));
title('CAN')

subplot(3,3,2)
normplot(mz_sq(:));
title({'Sqrt' , 'CMCC'})

subplot(3,3,3)
normplot(nz_sq(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_sq(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_sq(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_sq(:));
title('UK')

subplot(3,3,8)
normplot(oz_sq(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_sqrt.png'])

%% 4th-root
figure(7)
subplot(3,3,1)
normplot(cz_4th(:));
title('CAN')

subplot(3,3,2)
normplot(mz_4th(:));
title({'4th root' , 'CMCC'})

subplot(3,3,3)
normplot(nz_4th(:));
title('CNRM')

subplot(3,3,4)
normplot(gz_4th(:));
title('GFDL')

subplot(3,3,5)
normplot(iz_4th(:));
title('IPSL')

subplot(3,3,6)
normplot(uz_4th(:));
title('UK')

subplot(3,3,8)
normplot(oz_4th(:));
title('obsGLMM')
print('-dpng',[ppath 'Norm_all_clim_obsglm_4rt.png'])



