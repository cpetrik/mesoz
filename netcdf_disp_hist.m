% Read CMIP6 netcdfs
% GFDL-ESM4 Hist

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

%% phyc
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_phyc_onedeg_global_monthly_1850_2014.nc'])


%% diatoms zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_phydiat_onedeg_global_monthly_1850_2014.nc'])


%% diaz zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_phydiaz_onedeg_global_monthly_1850_2014.nc'])

%% zmicro z vint
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_zmicro-vint_onedeg_global_monthly_1850_2014.nc'])
