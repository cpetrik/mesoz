% CMIP6 output Historic 1951-2014
% 200m integrations of zmeso
% surface chl

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/hist/';
load([cpath 'can_hist_zmeso_200_monthly_1950_2014.mat'],'latitude','longitude',...
    'zmeso_200','mod_time','mod_yr');
load([cpath 'can_hist_surf_chl_monthly_1950_2014.mat'],'schl','time','yr','runs');

%% DIFFERENT LENGTHS OF TIME
%zmeso 1951-2014
%chl   1950-2014
time2 = time(runs);
[ti,id] = intersect(time2,mod_time);

%%
cc = schl(:,:,id);
cz = zmeso_200;
cz(cz(:)<0) = 0;
clat = latitude;
clon = longitude;

clear zmeso_200 schl latitude longitude

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/hist/';
load([npath 'cnrm_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon','zmeso_200');
load([npath 'cnrm_hist_surf_chl_monthly_1950_2014.mat'],'schl');

%%
nc = schl;
nz = zmeso_200;
nz(nz(:)<0) = 0;
nlat = lat;
nlon = lon;

clear zmeso_200 schl lat lon

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'ukesm_isimip_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon','zmeso_200');
load([upath 'ukesm_isimip_hist_surf_chl_monthly_1950_2014.mat'],'schl');

%%
uc = schl;
uz = zmeso_200;
uz(uz(:)<0) = 0;
ulat = lat;
ulon = lon;

clear zmeso_200 schl lat lon

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'ipsl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon','zmeso_200');
load([ipath 'ipsl_hist_surf_chl_monthly_1950_2014.mat'],'schl');

%%
ic = fliplr(schl);
iz = zmeso_200;
iz(iz(:)<0) = 0;
ilat = lat;
ilon = lon;

clear zmeso_200 schl lat lon

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
load([gpath 'gfdl_hist_zmeso_200_monthly_1950_2014.mat'],'lat','lon','zmeso_200');
load([gpath 'gfdl_hist_surf_chl_monthly_1950_2014.mat'],'schl');

%%
gc = schl;
gz = zmeso_200;
gz(gz(:)<0) = 0;
glat = lat;
glon = lon;

clear zmeso_200 schl lat lon

%% correlations
[ci,cj,ct] = size(cz);
vcc = reshape(cc,ci*cj,ct);
vcz = reshape(cz,ci*cj,ct);
cid = find(~isnan(vcz(:,1)));
vcc = vcc(cid,:);
vcz = vcz(cid,:);
ccorr = diag(corr(vcc',vcz'));
ccorr2 = NaN*ones(ci,cj);
ccorr2(cid) = ccorr;
%%
[ni,nj,nt] = size(nz);
vnc = reshape(nc,ni*nj,nt);
vnz = reshape(nz,ni*nj,nt);
nid = find(~isnan(vnz(:,1)));
vnc = vnc(nid,:);
vnz = vnz(nid,:);
ncorr = diag(corr(vnc',vnz'));
ncorr2 = NaN*ones(ni,nj);
ncorr2(nid) = ncorr;

[ui,uj,ut] = size(uz);
vuc = reshape(uc,ui*uj,ut);
vuz = reshape(uz,ui*uj,ut);
ucorr = diag(corr(vuc',vuz'));
ucorr = reshape(ucorr,ui,uj);

%%
[ii,ij,it] = size(iz);
vic = reshape(ic,ii*ij,it);
viz = reshape(iz,ii*ij,it);
icorr = diag(corr(vic',viz'));
icorr = reshape(icorr,ii,ij);

%%
[gi,gj,gt] = size(gz);
vgc = reshape(gc,gi*gj,gt);
vgz = reshape(gz,gi*gj,gt);
gcorr = diag(corr(vgc',vgz'));
gcorr = reshape(gcorr,gi,gj);
% test1 = xcorr(vgc,vgz);
% test2 = corr(vgc',vgz');
% test22 = diag(test2);
% test3 = corrcoef(vgc,vgz,'Rows','complete');

%%
save('corr_hist_zmeso_200_schl_monthly_1950_2014.mat',...
    'ccorr2','clat','clon','ncorr2','nlat','nlon','ucorr','ulat','ulon',...
    'icorr','ilat','ilon','gcorr2','glat','glon')

%% maps
clatlim=[-90 90];
clonlim=[-180 180];

[lat_u,lon_u] = meshgrid(ulat,ulon);
[lat_i,lon_i] = meshgrid(ilat,ilon);
[lat_g,lon_g] = meshgrid(glat,glon);

%% 
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(clat,clon,ccorr2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CAN','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(nlat,nlon,ncorr2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'CNRM','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_u,lon_u,ucorr)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'UKESM','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_i,lon_i,icorr)
cmocean('balance')
caxis([-1 1])
cb = colorbar('Position',[0.85 0.4 0.03 0.4],'orientation','vertical');
xlabel(cb,'corr mesoz-chl')
text(0.2,1.65,'IPSL','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_g,lon_g,gcorr2)
cmocean('balance')
caxis([-1 1])
text(0.2,1.65,'GFDL','HorizontalAlignment','center','FontWeight','bold')
load coast;
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% subplot('Position',[0.41 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(lat_g,lon_g,(zoo_200))
% cmocean('balance')
% caxis([-1 1])
% text(0.2,1.65,'COPEPOD mgC m^-^2','HorizontalAlignment','center','FontWeight','bold')
% load coast;
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_corr_hist_schl_zmeso_int200m_v2.png'])
 






