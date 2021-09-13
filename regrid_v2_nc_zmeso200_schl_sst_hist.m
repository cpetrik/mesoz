% CMIP6 output Historic 1951-2014
% 200m integrations of zmeso
% surface chl, SST
% take same time period 1951-2014 (most are 1950-2014)

clear all
close all

%% 1 degree
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[Lat1,Lon1] = meshgrid(lat1,lon1);

%% CAN ---------------------------------------------------------------
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso200_monthly_onedeg_1951_2014.mat'],...
    'zmeso200','time');

%%
yr_cz = time;

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
cz = zmeso200 * 12.01 * 1e3;
cz(cz(:)<0) = 0;

% [LAT,LON] = meshgrid(lat,lon);
% 
% figure
% pcolor(LAT)
% shading flat
% 
% figure
% pcolor(zmeso200(:,:,10))
% shading flat

%%
clear zmeso200 lat lon

%%
load([cpath 'can_hist_surf_chl_monthly_onedeg_1951_2014.mat'],'schl');
load([cpath 'can_hist_sst_monthly_onedeg_1951_2014.mat'],'sst',...
    'time','yr','runs');

%% 
yr_ct = yr(runs);

cchl = double(schl);
csst = double(sst);

%%
clear schl sst lat lon time yr runs 

%% CNRM ---------------------------------------------------------------
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso200_monthly_onedeg_1951_2014.mat'],...
    'zmeso200','yr','time');

%%
runs = 13:780;
yr_nz = time(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
nz = double(zmeso200(:,:,runs)) * 12.01 * 1e3;
nz(nz(:)<0) = 0;

%%
clear zmeso200 lat lon

%%
load([npath 'cnrm_hist_surf_chl_monthly_onedeg_1951_2014.mat'],'schl');
load([npath 'cnrm_hist_sst_monthly_onedeg_1951_2014.mat'],'sst',...
    'time','yr','runs');

%%
yr_nt = yr(runs);

nchl = double(schl);
nsst = double(sst);

%%
clear schl sst lat lon time yr runs 

%% UKESM ---------------------------------------------------------------
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'],...
    'lat','lon','zmeso_200','yr');

runs = (1980-768+1):1980;
yr_uz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
uz = double(zmeso_200(:,:,13:780)) * 12.01 * 1e3;
uz(uz(:)<0) = 0;

%%
clear zmeso_200 lat lon yr runs

%%
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([upath 'ukesm_isimip_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr');

runs = (1980-768+1):1980;
yr_ut = yr(runs);

uchl = double(schl(:,:,13:780));
usst = double(sst(:,:,13:780));

%%
clear schl sst time yr runs 

%% IPSL ---------------------------------------------------------------
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon',...
    'zmeso_200','yr');

runs = (1980-768+1):1980;
yr_iz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
iz = double(zmeso_200(:,:,13:780)) * 12.01 * 1e3;
iz(iz(:)<0) = 0;

%% this might need to be done because lon is whole numbers, not 0.5
% [~,~,it] = size(iz);
% [x,y] = size(Lat1);
% zI = NaN*ones(x,y,it);
% for t = 1:it
%     zI(:,:,t) = interp2(lat,lon,iz(:,:,t),Lat1,Lon1);
% end

%%
clear zmeso_200 lat lon yr runs

%%
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([ipath 'ipsl_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr','lat','lon');

runs = (1980-768+1):1980;
yr_it = yr(runs);

ichl = double(schl(:,:,13:780));
isst = double(sst(:,:,13:780));

%% this might need to be done because lon is whole numbers, not 0.5
% [ii,ij,it] = size(ichl);
% [x,y] = size(Lat1);
% cI = NaN*ones(x,y,it);
% tI = NaN*ones(x,y,it);
% for t = 1:it
%     cI(:,:,t) = interp2(lat,lon,ichl(:,:,t),Lat1,Lon1);
%     tI(:,:,t) = interp2(lat,lon,isst(:,:,t),Lat1,Lon1);
% end

%%
clear schl sst lat lon time yr runs 

%% GFDL ---------------------------------------------------------------
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'],...
    'zmeso_200','yr');

runs = (1980-768+1):1980;
yr_gz = yr(runs);

% Units
% meso zoo: from molC m-2 to mgC m-2
% 12.01 g C in 1 mol C
% 1e3 mg in 1 g
gz = double(zmeso_200(:,:,13:780)) * 12.01 * 1e3;
gz(gz(:)<0) = 0;

%%
clear zmeso_200 yr runs

%%
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],'schl');
load([gpath 'gfdl_hist_sst_monthly_1950_2014.mat'],'sst',...
    'yr');

runs = (1980-768+1):1980;
yr_gt = yr(runs);

gchl = double(schl(:,:,13:780));
gsst = double(sst(:,:,13:780));

%%
clear schl sst time yr runs

%% check orientations
figure(1)
pcolor(cz(:,:,1))
title('CAN zoo')

figure(2)
pcolor(cchl(:,:,4))
title('CAN chl')

figure(3)
pcolor(csst(:,:,46))
title('CAN sst')


figure(4)
pcolor(nz(:,:,11))
title('CNRM zoo')

figure(5)
pcolor(nchl(:,:,14))
title('CNRM chl')

figure(6)
pcolor(nsst(:,:,146))
title('CNRM sst')


figure(7)
pcolor(uz(:,:,21))
title('UK zoo')

figure(8)
pcolor(uchl(:,:,24))
title('UK chl')

figure(9)
pcolor(usst(:,:,246))
title('UK sst')


figure(10)
pcolor(iz(:,:,31))
title('IP zoo')

figure(11)
pcolor(ichl(:,:,34))
title('IP chl')

figure(12)
pcolor(isst(:,:,346))
title('IP sst')


figure(13)
pcolor(gz(:,:,41))
title('GF zoo')

figure(14)
pcolor(gchl(:,:,44))
title('GF chl')

figure(15)
pcolor(gsst(:,:,446))
title('GF sst')

%% flip as necessary
zC = fliplr(cz);
zG = gz;
zI = iz;
zN = fliplr(nz);
zU = fliplr(uz);

tC = fliplr(csst);
tG = gsst;
tI = isst;
tN = fliplr(nsst);
tU = usst;

cC = fliplr(cchl);
cG = gchl;
cI = fliplr(ichl);
cN = fliplr(nchl);
cU = fliplr(uchl);


%% save
save('/Volumes/MIP/Fish-MIP/CMIP6/Gridded_v2_hist_zmeso200_schl_sst_1951_2014.mat',...
    'cC','tC','zC','yr_cz','yr_ct',...
    'cN','tN','zN','yr_nz','yr_nt',...
    'cU','tU','zU','yr_uz','yr_ut',...
    'cI','tI','zI','yr_iz','yr_it',...
    'cG','tG','zG','yr_gz','yr_gt');




