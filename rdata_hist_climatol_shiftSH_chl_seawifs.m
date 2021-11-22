% Save obs chl used to force Stromberg model
% For calc linear regression

clear all
close all

%% Chl, SST, GLM zmeso, grid vars
opath ='/Volumes/MIP/Obs_Data/Chl/';
load([opath 'SeaWiFS_Mission_Climatology.mat'])

%% Reshape
[lat_c,lon_c] = meshgrid(lat,lon);

% climatologies
schl_all = chl_all;

%% check orientations
figure(1)
pcolor(schl_all); shading flat;
title('Obs')

figure(2)
pcolor(lat_c); shading flat;
title('Lat')

%% save same orientation
save('cmip6_hist_space_means_50yr_chl_sst_glmm100_same_orientation.mat',...
    'schl_all','-append');

%% Vectorize, put in ABC order
% convert to same units
%models kg/m3 -> mg/m3 CAN, GFDL, UK
%models g/m3 -> mg/m3  CNRM, IPSL
%obsglm mg/m3 ?
mod_chl = (schl_all(:));

%% all clim
lat = lat_c(:);
lon = lon_c(:);

matc(:,1) = lat;
matc(:,2) = lon;
matc(:,3) = mod_chl;

tchl = array2table(matc,'VariableNames',...
    {'Lat','Lon','seawifs'});
writetable(tchl,'hist_chl_seawifs_all_clim.csv')

%%
matc_seawifs = matc;
save('hist_chl_sst_model_obsglm100_climatols.mat','matc_seawifs','-append')








