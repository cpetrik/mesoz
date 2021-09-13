clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp585/';
path2='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/ssp585/';
path3='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
path4='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';
path5='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
%%
ncdisp([fpath 'zmeso200_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_2015_2100_onedeg.nc'])
%%
ncdisp([fpath 'tos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'])
%%
ncdisp([fpath 'chlos_Omon_CanESM5-CanOE_ssp585_r1i1p2f1_gn_201501-210012_onedeg.nc'])
%%
ncdisp([path2 'zmeso200_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_2015_2100_onedeg.nc'])
%%
ncdisp([path2 'tos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012_onedeg.nc'])
%%
ncdisp([path2 'chlos_Omon_CNRM-ESM2-1_ssp585_r1i1p1f2_gn_201501-210012_onedeg.nc'])





