% Read CMIP6 netcdfs
% GFDL-ESM4 Hist
% Epc100
% Regular grid from ESGF

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/GFDL_CMIP6/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/Phase3/quarter_degree_vetting/';

%% expc
ncdisp([fpath 'epc100_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc'])

%%
% Dimensions:
% time = 240   (UNLIMITED)
% lat  = 180
% lon  = 360
% bnds = 2

% Variables:
% epc100
% Size:       360x180x240
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
long_name     = 'Downward Flux of Particulate Organic Carbon';
units         = 'mol m-2 s-1';
missing_value = 1.000000020040877e+20;
FillValue    = 1.000000020040877e+20;
% cell_methods  = 'area: mean where sea time: mean'
% cell_measures = 'area: areacello'
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
% interp_method = 'conserve_order1'
% original_name = 'epc100'
% comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'

%%
ncid = netcdf.open([fpath 'epc100_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

epc100(epc100>1e19) = nan;

%% Time
mepc100 = double(mean(epc100,3));

%%
save([fpath 'epc_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.mat'],...
    'epc','time',...
    'long_name','units','lat','lon');

%% Convert units
%molC m-2 s-1 to kgC m-2
%12.0gC in 1 molC, 1e-3kg per g, mult by seconds in a day? month?
mepc = mepc100 * 12.0 * 1e-3 * 3600 * 24 * 30;

%% Map
[gLAT,gLON] = meshgrid(lat,lon);

clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%% MAPS
%1 - expc
figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,mepc)
%cmocean('matter')
colormap('jet')
caxis([6e-5 0.007])
colorbar
title('epc-100')
%text(-2.5,2.25,num2str(qint_poc),'HorizontalAlignment','left')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_GFDL_CMIP6_1deg_from_ESGF_jet_epc100.png'])




