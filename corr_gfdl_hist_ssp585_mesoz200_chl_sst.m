% Save correlations w/mesoz Hist and SSP 
% GFDL mesoz prepared by ISIMIP team

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';
ppath='/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%%
load([hpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_sst_monthly_1950_2014.mat'])
load([hpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'])

%%
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

% Hist time
Hyr = yr(runs);
Hzm = double(zmeso_200(:,:,(end-599):end));
Hsst = double(sst(:,:,(end-599):end));
Hschl = double(schl(:,:,(end-599):end));

[lat_g,lon_g] = meshgrid(lat,lon);

clear zmeso_200 sst schl

%% check orientation
figure
pcolor(squeeze(Hzm(:,:,50)))
shading flat
title('z')

figure
pcolor(squeeze(Hschl(:,:,50)))
shading flat
title('chl')

figure
pcolor(squeeze(Hsst(:,:,50)))
shading flat
title('sst')

figure
pcolor(lat_g)
shading flat
title('lat')

%% SSP 
load([spath 'gfdl_ssp585_zmeso_200_monthly_2015_2100.mat'])
load([spath 'gfdl_ssp585_sst_monthly_2015_2100.mat'])
load([spath 'gfdl_ssp585_surf_chl_monthly_2015_2100.mat'])

%% 
zmeso_200(zmeso_200>=1e19) = NaN;
sst(sst>=1e19) = NaN;
schl(schl>=1e19) = NaN;

Fyr = yr;
Fzm = double(zmeso_200(:,:,(end-599):end));
Fsst = double(sst(:,:,(end-599):end));
Fschl = double(schl(:,:,(end-599):end));

[lat_g,lon_g] = meshgrid(lat,lon);

clear zmeso_200 sst schl

%% check orientation
figure
pcolor(squeeze(Fzm(:,:,50)))
shading flat
title('z')

figure
pcolor(squeeze(Fschl(:,:,50)))
shading flat
title('chl')

figure
pcolor(squeeze(Fsst(:,:,50)))
shading flat
title('sst')

figure
pcolor(lat_g)
shading flat
title('lat')

%% Hist correlations by grid cell
[ni,nj,nt] = size(Hzm);
vc = reshape(Hschl,ni*nj,nt);
vz = reshape(Hzm,ni*nj,nt);
vt = reshape(Hsst,ni*nj,nt);

%remove land cells
id = find(~isnan(vz(:,1)));
vc = vc(id,:);
vz = vz(id,:);
vt = vt(id,:);

ccorr = diag(corr(vc',vz'));
tcorr = diag(corr(vt',vz'));

%gccorr = reshape(ccorr,ni,nj);
%gtcorr = reshape(tcorr,ni,nj);
gccorr = nan(ni,nj);
gtcorr = nan(ni,nj);

gccorr(id) = ccorr;
gtcorr(id) = tcorr;

%% SSP correlations by grid cell
clear vc vz vt ni nj nt id ccorr tcorr

[ni,nj,nt] = size(Fzm);
vc = reshape(Fschl,ni*nj,nt);
vz = reshape(Fzm,ni*nj,nt);
vt = reshape(Fsst,ni*nj,nt);

%remove land cells
id = find(~isnan(vz(:,1)));
vc = vc(id,:);
vz = vz(id,:);
vt = vt(id,:);

ccorr2 = diag(corr(vc',vz'));
tcorr2 = diag(corr(vt',vz'));

%gccorr = reshape(ccorr,ni,nj);
%gtcorr = reshape(tcorr,ni,nj);
gccorr2 = nan(ni,nj);
gtcorr2 = nan(ni,nj);

gccorr2(id) = ccorr2;
gtcorr2(id) = tcorr2;

%% Plots
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

figure(1)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,gccorr)
cmocean('balance')
caxis([-1 1])
title('chl')
text(0.2,1.65,'GFDL Hist','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,gtcorr)
cmocean('balance')
caxis([-1 1])
title('sst')
text(0.2,1.65,'GFDL Hist','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,gccorr2)
cmocean('balance')
caxis([-1 1])
title('chl')
text(0.2,1.65,'GFDL SSP585','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,gtcorr2)
cmocean('balance')
caxis([-1 1])
title('sst')
text(0.2,1.65,'GFDL SSP585','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%
save([fpath 'gfdl_hist_ssp585_corrs_grid_zmeso200_schl_sst.mat'],'lat_g','lon_g',...
    'gccorr','gtcorr','gccorr2','gtcorr2');
