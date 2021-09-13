% Calculate different skill metrics for each ESM
% log transform biomass
% shift southern hemisphere by 6 mo (Summer = DJF)

clear all
close all

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all');

chl_all(chl_all(:)<0) = 0;

cchl_all = chl_all;
csst_all = sst_all;

clear chl_all sst_all

%% CMCC
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/CMCC/hist/';
load([mpath 'cmcc_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all');

chl_all(chl_all(:)<0) = 0;

mchl_all = chl_all;
msst_all = sst_all;

clear chl_all sst_all

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all');

chl_all(chl_all(:)<0) = 0;

nchl_all = chl_all;
nsst_all = sst_all;

clear chl_all sst_all

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all');

chl_all(chl_all(:)<0) = 0;

uchl_all = chl_all;
usst_all = sst_all;

clear chl_all sst_all

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all');

chl_all(chl_all(:)<0) = 0;

ichl_all = chl_all;
isst_all = sst_all;

clear chl_all sst_all

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_chl_sst_onedeg_climatol_1965_2014.mat'],...
    'chl_all','sst_all','lat','lon');

chl_all(chl_all(:)<0) = 0;

gchl_all = chl_all;
gsst_all = sst_all;

clear chl_all sst_all

%% Chl, SST, GLM zmeso, grid vars
opath ='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/glm_Ryan/';
load([opath 'glm_obs_chl.mat']);
load([opath 'glm_obs_sst.mat']);
load([opath 'glm_obs_grid.mat'])

% Reshape
[ni,nj] = size(gchl_all);
glmc_mo = reshape(glmobschl,ni,nj,12);
glmt_mo = reshape(glmobssst,ni,nj,12);
lat_g = reshape(Lat,ni,nj);
lon_g = reshape(Lon,ni,nj);

% climatologies
ochl_all = nanmean(glmc_mo,3);
osst_all = nanmean(glmt_mo,3);

%% check orientations
figure(1)
pcolor(cchl_all); shading flat;
title('CAN')

figure(7)
pcolor(mchl_all); shading flat;
title('CMCC')

figure(2)
pcolor(nchl_all); shading flat;
title('CNRM')

figure(3)
pcolor(ichl_all); shading flat;
title('IPSL')

figure(4)
pcolor(gchl_all); shading flat;
title('GFDL')

figure(5)
pcolor(uchl_all); shading flat;
title('UK')

figure(6)
pcolor(ochl_all); shading flat;
title('Obs')

%% check orientations sst
close all
figure(1)
pcolor(csst_all); shading flat;
title('CAN')

figure(7)
pcolor(msst_all); shading flat;
title('CMCC')

figure(2)
pcolor(nsst_all); shading flat;
title('CNRM')

figure(3)
pcolor(isst_all); shading flat;
title('IPSL')

figure(4)
pcolor(gsst_all); shading flat;
title('GFDL')

figure(5)
pcolor(usst_all); shading flat;
title('UK')

figure(6)
pcolor(osst_all); shading flat;
title('Obs')

%% Flip
%UK, CNRM, CAN
cchl_all = fliplr(cchl_all);
mchl_all = fliplr(mchl_all);
nchl_all = fliplr(nchl_all);
ichl_all = fliplr(ichl_all);
uchl_all = fliplr(uchl_all);

csst_all = fliplr(csst_all);
msst_all = fliplr(msst_all);
nsst_all = fliplr(nsst_all);

%% Vectorize, put in ABC order
% convert to same units
%models kg/m3 -> mg/m3 CAN, GFDL, UK
%models g/m3 -> mg/m3  CNRM, IPSL
%obsglm mg/m3 
mod_chl(:,1) = (cchl_all(:)) * 1e6;
mod_chl(:,2) = (mchl_all(:)) * 1e6;
mod_chl(:,3) = (nchl_all(:)) * 1e3;
mod_chl(:,4) = (gchl_all(:)) * 1e6;
mod_chl(:,5) = (ichl_all(:)) * 1e3;
mod_chl(:,6) = (uchl_all(:)) * 1e6;
mod_chl(:,7) = (ochl_all(:));

mod_sst(:,1) = (csst_all(:));
mod_sst(:,2) = (msst_all(:));
mod_sst(:,3) = (nsst_all(:));
mod_sst(:,4) = (gsst_all(:));
mod_sst(:,5) = (isst_all(:));
mod_sst(:,6) = (usst_all(:));
mod_sst(:,7) = (osst_all(:));

%% all clim
lat = lat_g(:);
lon = lon_g(:);

matc(:,1) = lat;
matc(:,2) = lon;
matc(:,3:9) = mod_chl;

matt(:,1) = lat;
matt(:,2) = lon;
matt(:,3:9) = mod_sst;

tchl = array2table(matc,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(tchl,'hist_chl_model_obsglm_all_clim.csv')

tsst = array2table(matt,'VariableNames',...
    {'Lat','Lon','CAN','CMCC','CNRM','GFDL','IPSL','UK','obsGLM'});
writetable(tsst,'hist_sst_model_obsglm_all_clim.csv')

%%
save('hist_chl_sst_model_obsglm_climatols.mat','matc','matt')








