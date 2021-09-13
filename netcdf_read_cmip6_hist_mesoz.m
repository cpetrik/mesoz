% Check units of CMIP6 Historic netcdfs

clear all
close all

%% CAN Meso Zoop zall
fpath='/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
ncdisp([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_201101-201412.nc'])

ncid = netcdf.open([fpath 'zmeso_Omon_CanESM5-CanOE_historical_r1i1p2f1_gn_201101-201412.nc'],'NC_NOWRITE');
    
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

zmeso(zmeso >= 1.00e+20) = NaN;
%Surf at Jan 2011
zmeso_can = double(squeeze(zmeso(:,:,1,1)));
clear zmeso

cLAT = latitude;
cLON = longitude;

%% CNRM
cpath='/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
ncdisp([cpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'])
                         
%
ncid = netcdf.open([cpath 'zmeso_Omon_CNRM-ESM2-1_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of zmeso so smaller memory
yr2 = ((time+1)/365)+1850;
runs1 = find(yr2>2011);
z1 = 1:5;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 294 length(z1) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end

%% Surf at Jan 2011
zmeso_cnrm = double(squeeze(zmeso(:,:,1,1)));
clear zmeso

nLAT = lat;
nLON = lon;

%% GFDL
gpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
ncdisp([gpath 'zmeso_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc'])
                         
%% 
ncid = netcdf.open([gpath 'zmeso_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
zmeso(zmeso >= 1.00e+20) = NaN;

%% Surf at Jan 
zmeso_gfdl = double(squeeze(zmeso(:,:,1,1)));
clear zmeso

[gLAT,gLON] = meshgrid(lat,lon);

%% IPSL
ipath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
ncdisp([ipath 'zmeso_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_195001-201412.nc'])
ncid = netcdf.open([ipath 'zmeso_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_195001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of zmeso so smaller memory
yr2 = ((time+1)/365)+1850;
runs1 = find(yr2>2011);
z1 = 1:5;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[362 332 length(z1) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end

%% Surf at Jan 2011
zmeso_ipsl = double(squeeze(zmeso(:,:,1,1)));
clear zmeso

iLAT = double(nav_lat);
iLON = double(nav_lon);

%% UK
upath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
ncdisp([upath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'])
                         
%
ncid = netcdf.open([upath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of zmeso so smaller memory
yr1 = ((time+1)/360)+1850;
runs1 = find(yr1>2011);
z1 = 1:5;
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs1(1)-1],[360 330 length(z1) length(runs1)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end

%% Surf at Jan 2011
zmeso_uk = double(squeeze(zmeso(:,:,1,1)));
clear zmeso
zmeso_uk(zmeso_uk<0) = 0;

uLAT = double(latitude);
uLON = double(longitude);

%%
edg = [-20 -15 -10 -5 -4 -3 -2 -1 0 1];

figure(2)
subplot(3,2,1)
histogram(log10(zmeso_cnrm(:)),edg)
title('CNRM')
subplot(3,2,2)
histogram(log10(zmeso_can(:)),edg)
title('CAN')
subplot(3,2,3)
histogram(log10(zmeso_gfdl(:)),edg)
title('GFDL')
subplot(3,2,4)
histogram(log10(zmeso_ipsl(:)),edg)
title('IPSL')
subplot(3,2,5)
histogram(log10(zmeso_uk(:)),edg)
title('UK')
print('-dpng','Bar_all_hist_clim_molCm3_Jan2011.png')

%%
clatlim=[-90 90];
clonlim=[-180 180];

%% Maps
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(cLAT,cLON,zmeso_can)
cmocean('tempo')
colorbar
caxis([0 3e-3])
title('CAN Jan 2011')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(nLAT,nLON,zmeso_cnrm)
cmocean('tempo')
colorbar
caxis([0 3e-3])
%cb = colorbar('Position',[0.3 0.475 0.4 0.03],'orientation','horizontal');
%xlabel(cb,'log_1_0 mesoz (mg m^-^3)')
title('CNRM Jan 2011')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,zmeso_gfdl)
cmocean('tempo')
colorbar
caxis([0 3e-3])
title('GFDL Jan 2011')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(iLAT,iLON,zmeso_ipsl)
cmocean('tempo')
colorbar
caxis([0 3e-3])
title('IPSL Jan 2011')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng','Map_all_hist_clim_molCm3_Jan2011.png')

%%
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(uLAT,uLON,zmeso_uk)
cmocean('tempo')
colorbar
caxis([0 3e-3])
title('UK Jan 2011')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
