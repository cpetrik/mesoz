% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% MLD on regular grid

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/mlotst/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% mld
ncdisp([fpath 'mlotst_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc'])

%%
long_name     = 'Ocean Mixed Layer Thickness Defined by Sigma T';
units         = 'm';
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
standard_name = 'ocean_mixed_layer_thickness_defined_by_sigma_t';
% interp_method = 'conserve_order1'
% original_name = 'mlotst'
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'
% time units    = 'days since 1850-01-01 00:00:00'
% calendar      = 'noleap'

%% time slices
ttex = {'201501-203412', '203501-205412', '205501-207412','207501-209412',...
    '209501-210012'};

%%
mlot =[];
tall = [];
for t = 1:length(ttex)
ncid = netcdf.open([fpath 'mlotst_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_',ttex{t},'.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

mlot = cat(3,mlot,mlotst);
tall = [tall; time];

end

%% Time
yr = ((tall+1)/365)+1850;
runs = find(yr>2051);% & yr<=2100);
mld = mlot(:,:,runs);

%%
save([spath 'gfdl_ssp585_mld_monthly_2051_2100.mat'],'mld','yr',...
    'long_name','standard_name','units','lat','lon',...
    'runs');





