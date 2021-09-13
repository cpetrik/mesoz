% CMIP6 output Historic 1951-2014
% 200m integrations of zmeso
% surface chl, SST

clear all
close all

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[Lat1,Lon1] = meshgrid(lat1,lon1);

%% CAN ---------------------------------------------------------------
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso_200_monthly_1950_2014.mat'],'latitude','longitude',...
    'zmeso_200','mod_yr');

yr_cz = mod_yr;

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
cz = double(zmeso_200) * 12.01 * 1e3;
cz(cz(:)<0) = 0;

% longitude is shifted 0-360 and grid is not regulat (tripolar?)
max(longitude(:))
min(longitude(:))

test = longitude; %-360;
id=find(test > 180);
test(id)=test(id)-360;
LON = double(test);
LAT = latitude;

%
[ci,cj,ct] = size(cz);
[x,y] = size(Lat1);
zC = NaN*ones(x,y,ct);
for t = 1:ct
    zC(:,:,t) = griddata(LAT,LON,cz(:,:,t),Lat1,Lon1);
end

%
clear zmeso_200 latitude longitude

%
load([cpath 'can_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([cpath 'can_hist_sst_monthly_1950_2014.mat'],'sst',...
    'time','yr','runs','latitude','longitude');

% DIFFERENT LENGTHS OF TIME
%zmeso 1951-2014
%chl   1950-2014
yr_ct = yr(runs);
id = find(yr_ct > 1951);

schl(schl(:)>=1e19) = NaN;
sst(sst(:)>=1e19) = NaN;

cchl = double(schl(:,:,id));
csst = double(sst(:,:,id));

% longitude is shifted 0-360 and grid is not regulat (tripolar?)
max(longitude(:))
min(longitude(:))

test = longitude; %-360;
id=find(test > 180);
test(id)=test(id)-360;
LON = double(test);
LAT = latitude;

%
[ci,cj,ct] = size(cchl);
[x,y] = size(Lat1);
cC = NaN*ones(x,y,ct);
tC = NaN*ones(x,y,ct);
for t = 1:ct
    cC(:,:,t) = griddata(LAT,LON,cchl(:,:,t),Lat1,Lon1);
    tC(:,:,t) = griddata(LAT,LON,csst(:,:,t),Lat1,Lon1);
end

%
clear schl sst latitude longitude time yr runs LAT LON

% CNRM ---------------------------------------------------------------
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon',...
    'zmeso_200','yr1','yr2');

yr_nz = [yr1; yr2];

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
nz = double(zmeso_200) * 12.01 * 1e3;
nz(nz(:)<0) = 0;

%
[ni,nj,nt] = size(nz);
[x,y] = size(Lat1);
zN = NaN*ones(x,y,nt);
for t = 1:nt
    zN(:,:,t) = griddata(lat,lon,nz(:,:,t),Lat1,Lon1);
end

%
clear zmeso_200 lat lon

%%
load([npath 'cnrm_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([npath 'cnrm_hist_sst_monthly_1950_2014.mat'],'sst',...
    'time','yr','runs','lat','lon');

yr_nt = yr(runs);

schl(schl(:)>=1e19) = NaN;
sst(sst(:)>=1e19) = NaN;

nchl = double(schl);
nsst = double(sst);

%%
[ni,nj,nt] = size(nchl);
[x,y] = size(Lat1);
cN = NaN*ones(x,y,nt);
tN = NaN*ones(x,y,nt);
for t = 1:nt
    cN(:,:,t) = griddata(lat,lon,nchl(:,:,t),Lat1,Lon1);
    tN(:,:,t) = griddata(lat,lon,nsst(:,:,t),Lat1,Lon1);
end

%
clear schl sst lat lon time yr runs 


%% UKESM ---------------------------------------------------------------
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'],...
    'lat','lon','zmeso_200','yr','runs');

yr_uz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
uz = double(zmeso_200) * 12.01 * 1e3;
uz(uz(:)<0) = 0;

%%
[ui,uj,ut] = size(uz);
[x,y] = size(Lat1);
zU = NaN*ones(x,y,ut);
for t = 1:ut
    zU(:,:,t) = interp2(lat,lon,uz(:,:,t),Lat1,Lon1);
end

%%
clear zmeso_200 lat lon yr runs

%%
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([upath 'ukesm_isimip_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr','runs','lat','lon');

yr_ut = yr(runs);

uchl = double(schl);
usst = double(sst);

%%
[ui,uj,ut] = size(uchl);
[x,y] = size(Lat1);
cU = NaN*ones(x,y,ut);
tU = NaN*ones(x,y,ut);
for t = 1:ut
    cU(:,:,t) = interp2(lat,lon,uchl(:,:,t),Lat1,Lon1);
    tU(:,:,t) = interp2(lat,lon,usst(:,:,t),Lat1,Lon1);
end

%%
clear schl sst lat lon time yr runs 

%% IPSL ---------------------------------------------------------------
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon',...
    'zmeso_200','yr','runs');

yr_iz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
iz = double(zmeso_200) * 12.01 * 1e3;
iz(iz(:)<0) = 0;

%%
[ii,ij,it] = size(iz);
[x,y] = size(Lat1);
zI = NaN*ones(x,y,it);
for t = 1:it
    zI(:,:,t) = interp2(lat,lon,iz(:,:,t),Lat1,Lon1);
end

%%
clear zmeso_200 lat lon yr runs

%%
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([ipath 'ipsl_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr','runs','lat','lon');

yr_it = yr(runs);

ichl = double(schl);
isst = double(sst);

%%
[ii,ij,it] = size(ichl);
[x,y] = size(Lat1);
cI = NaN*ones(x,y,it);
tI = NaN*ones(x,y,it);
for t = 1:it
    cI(:,:,t) = interp2(lat,lon,ichl(:,:,t),Lat1,Lon1);
    tI(:,:,t) = interp2(lat,lon,isst(:,:,t),Lat1,Lon1);
end

%%
clear schl sst lat lon time yr runs 

%% GFDL ---------------------------------------------------------------
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon',...
    'zmeso_200','yr','runs');

yr_gz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
gz = double(zmeso_200) * 12.01 * 1e3;
gz(gz(:)<0) = 0;

%%
[gi,gj,gt] = size(gz);
[x,y] = size(Lat1);
zG = NaN*ones(x,y,gt);
for t = 1:gt
    zG(:,:,t) = interp2(lat,lon,gz(:,:,t),Lat1,Lon1);
end

%%
clear zmeso_200 lat lon yr runs

%%
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([gpath 'gfdl_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr','runs','lat','lon');

yr_gt = yr(runs);

gchl = double(schl);
gsst = double(sst);

%%
[gi,gj,gt] = size(gchl);
[x,y] = size(Lat1);
cG = NaN*ones(x,y,gt);
tG = NaN*ones(x,y,gt);
for t = 1:gt
    cG(:,:,t) = interp2(lat,lon,gchl(:,:,t),Lat1,Lon1);
    tG(:,:,t) = interp2(lat,lon,gsst(:,:,t),Lat1,Lon1);
end

%%
clear schl sst lat lon time yr runs

%% check orientations
figure(1)
pcolor(zC(:,:,1))
title('CAN zoo')

figure(2)
pcolor(cC(:,:,4))
title('CAN chl')

figure(3)
pcolor(tC(:,:,46))
title('CAN sst')


figure(4)
pcolor(zN(:,:,11))
title('CNRM zoo')

figure(5)
pcolor(cN(:,:,14))
title('CNRM chl')

figure(6)
pcolor(tN(:,:,146))
title('CNRM sst')


figure(7)
pcolor(zU(:,:,21))
title('UK zoo')

figure(8)
pcolor(cU(:,:,24))
title('UK chl')

figure(9)
pcolor(tU(:,:,246))
title('UK sst')


figure(10)
pcolor(zI(:,:,31))
title('IP zoo')

figure(11)
pcolor(cI(:,:,34))
title('IP chl')

figure(12)
pcolor(tI(:,:,346))
title('IP sst')


figure(13)
pcolor(zG(:,:,41))
title('GF zoo')

figure(14)
pcolor(cG(:,:,44))
title('GF chl')

figure(15)
pcolor(tG(:,:,446))
title('GF sst')

%% flip as necessary
tU = fliplr(tU);
cI = fliplr(cI);

%% take same time period 1951-2014
% most are 1950-2014
tid = 13:780;
zG = zG(:,:,tid);
zI = zI(:,:,tid);
zN = zN(:,:,tid);
zU = zU(:,:,tid);

tG = tG(:,:,tid);
tI = tI(:,:,tid);
tN = tN(:,:,tid);
tU = tU(:,:,tid);

cG = cG(:,:,tid);
cI = cI(:,:,tid);
cN = cN(:,:,tid);
cU = cU(:,:,tid);

%%
yr_ct = yr_ct(tid);
yr_gt = yr_gt(tid);
yr_it = yr_it(tid);
yr_nt = yr_nt(tid);
yr_ut = yr_ut(tid);

yr_gz = yr_gz(tid);
yr_iz = yr_iz(tid);
yr_nz = yr_nz(tid);
yr_uz = yr_uz(tid);


%% save
save('/Volumes/MIP/Fish-MIP/CMIP6/Gridded_hist_zmeso200_schl_sst_1951_2014.mat',...
    'cC','tC','zC','yr_cz','yr_ct',...
    'cN','tN','zN','yr_nz','yr_nt',...
    'cU','tU','zU','yr_uz','yr_ut',...
    'cI','tI','zI','yr_iz','yr_it',...
    'cG','tG','zG','yr_gz','yr_gt');




