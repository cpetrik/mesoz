% Read CMIP6 netcdfs
% IPSL Hist
% Vert int top 200 m

clear all
close all

%fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

%% Meso Zoop zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%% 
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);
z200 = find(olevel <= 200);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1], [360 180 length(z200) length(runs)]);
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thkcello_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% check if either are L-R flipped (lat,lon) 
whos thkcello zmeso
ttest = squeeze(thkcello(:,:,1,100));
ztest = squeeze(zmeso(:,:,1,100));

figure
pcolor(ztest)
shading flat

figure
pcolor(ttest)
shading flat

%% Integrate top 200 m
zmeso_200 = squeeze(nansum((zmeso.*thkcello),3));

%%
whos zmeso_200
ztest2 = squeeze(zmeso_200(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'],'zmeso_200','yr',...
    'long_name','standard_name','units_orig','units_vint','lat','lon',...
    'runs','z200','olevel');





