
% Read CMIP6 netcdfs
% CanOE SSP 126
% Vert int top 200 m

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/ssp126/';

%% Meso Zoop zall
ncdisp([fpath 'zmeso_Omon_CanESM5-CanOE_ssp126_r1i1p2f1_gn_201501-210012.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x291x45x1032
%Dimensions: i,j,time
%time units = 'days since 1850-01-01 '
                         
%% 
ncid = netcdf.open([fpath 'zmeso_Omon_CanESM5-CanOE_ssp126_r1i1p2f1_gn_201501-210012.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time)/365)+1850;
z200 = find(lev <= 200);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,0], [360 291 length(z200) length(time)]);
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
ncdisp([fpath 'thkcello_Ofx_CanESM5-CanOE_ssp126_r1i1p2f1_gn.nc'])

%%
tcid = netcdf.open([fpath 'thkcello_Ofx_CanESM5-CanOE_ssp126_r1i1p2f1_gn.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0],[360 291 length(z200)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

thkcello = repmat(thkcello,1,1,1,length(time));

%% Make sure oriented the same
whos zmeso thkcello
ztest = squeeze(zmeso(:,:,1,150));
ttest = squeeze(thkcello(:,:,1,150));

figure
pcolor(ztest)
figure
pcolor(ttest)

%% Integrate top 200 m
zmeso_200 = squeeze(nansum((zmeso.*thkcello),3));

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'can_ssp126_zmeso_200_monthly_2015_2100.mat'],'zmeso_200','yr',...
    'long_name','standard_name','units_orig','units_vint','latitude','longitude',...
    'time','z200','lev');





