% CMIP6 output 
% mesoz 200m integrations
% Map ESM future diffs 2090s vs. 1990s

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

%% CAN
cpath = '/Volumes/MIP/Fish-MIP/CMIP6/CanESM5/';
load([cpath 'can_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% CNRM
npath = '/Volumes/MIP/Fish-MIP/CMIP6/CNRM-ESM2-1/';
load([npath 'cnrm_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% UKESM
upath = '/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/';
load([upath 'ukesm_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

%% IPSL
ipath = '/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
load([ipath 'ipsl_hist_ssp585_space_means_zmeso200_schl_sst.mat']); 

%% GFDL
gpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';
load([gpath 'gfdl_hist_ssp585_space_means_zmeso200_schl_sst.mat']);

load([gpath 'hist/gfdl_hist_zmeso200_onedeg_climatol_1965_2014.mat'],...
    'lat','lon');
[lat_g,lon_g] = meshgrid(lat,lon);

%% Convert all zoo units to mgC/m2
%all models in molC: 12.01 g C in 1 mol C
%1e3 mg in 1 g
HCzm = HCzm90 * 12.01 * 1e3;
HNzm = HNzm90 * 12.01 * 1e3;
HGzm = HGzm90 * 12.01 * 1e3;
HIzm = HIzm90 * 12.01 * 1e3;
HUzm = HUzm90 * 12.01 * 1e3;

FCzm = FCzm90 * 12.01 * 1e3;
FNzm = FNzm90 * 12.01 * 1e3;
FGzm = FGzm90 * 12.01 * 1e3;
FIzm = FIzm90 * 12.01 * 1e3;
FUzm = FUzm90 * 12.01 * 1e3;

%chl in kg m-3', put in g m-3 
%chl in g m-3', put in mg m-3 (CNRM & IPSL)
HCchl = HCchl90 * 1e6;
HGchl = HGchl90 * 1e6;
HUchl = HUchl90 * 1e6;
HNchl = HNchl90 * 1e3;
HIchl = HIchl90 * 1e3;

FCchl = FCchl90 * 1e6;
FGchl = FGchl90 * 1e6;
FUchl = FUchl90 * 1e6;
FNchl = FNchl90 * 1e3;
FIchl = FIchl90 * 1e3;

%% hist z
figure
pcolor(HCzm)
shading flat
title('CAN')

figure
pcolor(HNzm)
shading flat
title('CNRM')

figure
pcolor(HIzm)
shading flat
title('IP')

figure
pcolor(HGzm)
shading flat
title('GF')

figure
pcolor(HUzm)
shading flat
title('UK')

%% ssp z
close all
figure
pcolor(FCzm)
shading flat
title('CAN')

figure
pcolor(FNzm)
shading flat
title('CNRM')

figure
pcolor(FIzm)
shading flat
title('IP')

figure
pcolor(FGzm)
shading flat
title('GF')

figure
pcolor(FUzm)
shading flat
title('UK')

%% hist chl
close all
figure
pcolor(HCchl)
shading flat
title('CAN')

figure
pcolor(HNchl)
shading flat
title('CNRM')

figure
pcolor(HIchl)
shading flat
title('IP')

figure
pcolor(HGchl)
shading flat
title('GF')

figure
pcolor(HUchl)
shading flat
title('UK')

%% ssp chl
close all
figure
pcolor(FCchl)
shading flat
title('CAN')

figure
pcolor(FNchl)
shading flat
title('CNRM')

figure
pcolor(FIchl)
shading flat
title('IP')

figure
pcolor(FGchl)
shading flat
title('GF')

figure
pcolor(FUchl)
shading flat
title('UK')

%% hist sst
close all
figure
pcolor(HCsst90)
shading flat
title('CAN')

figure
pcolor(HNsst90)
shading flat
title('CNRM')

figure
pcolor(HIsst90)
shading flat
title('IP')

figure
pcolor(HGsst90)
shading flat
title('GF')

figure
pcolor(HUsst90)
shading flat
title('UK')

%% ssp sst90
close all
figure
pcolor(FCsst90)
shading flat
title('CAN')

figure
pcolor(FNsst90)
shading flat
title('CNRM')

figure
pcolor(FIsst90)
shading flat
title('IP')

figure
pcolor(FGsst90)
shading flat
title('GF')

figure
pcolor(FUsst90)
shading flat
title('UK')

%% lat
close all
figure
pcolor(lat_g)
shading flat
title('lat')

%% flip as needed
close all
HGzm = fliplr(HGzm);
FGzm = fliplr(FGzm);
FUzm = fliplr(FUzm);
HGchl = fliplr(HGchl);
FGchl = fliplr(FGchl);
HIsst90 = fliplr(HIsst90);
HGsst90 = fliplr(HGsst90);
HUsst90 = fliplr(HUsst90);
FIsst90 = fliplr(FIsst90);
FGsst90 = fliplr(FGsst90);
FUsst90 = fliplr(FUsst90);
lat_g = fliplr(lat_g);

%% Diffs 
diff_ct = FCsst90 - HCsst90;
diff_nt = FNsst90 - HNsst90;
diff_gt = FGsst90 - HGsst90;
diff_it = FIsst90 - HIsst90;
diff_ut = FUsst90 - HUsst90;

%% Pdiffs
pdiff_cc = (FCchl - HCchl) ./ HCchl;
pdiff_nc = (FNchl - HNchl) ./ HNchl;
pdiff_gc = (FGchl - HGchl) ./ HGchl;
pdiff_ic = (FIchl - HIchl) ./ HIchl;
pdiff_uc = (FUchl - HUchl) ./ HUchl;

pdiff_cz = (FCzm - HCzm) ./ HCzm;
pdiff_nz = (FNzm - HNzm) ./ HNzm;
pdiff_gz = (FGzm - HGzm) ./ HGzm;
pdiff_iz = (FIzm - HIzm) ./ HIzm;
pdiff_uz = (FUzm - HUzm) ./ HUzm;

%% Fix seam
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);
[~,~,diff_nt2] = cyclic_map_seam(lat_g,lon_g,diff_nt);
[~,~,diff_gt2] = cyclic_map_seam(lat_g,lon_g,diff_gt);
[~,~,diff_it2] = cyclic_map_seam(lat_g,lon_g,diff_it);
[~,~,diff_ut2] = cyclic_map_seam(lat_g,lon_g,diff_ut);

[~,~,pdiff_cc2] = cyclic_map_seam(lat_g,lon_g,pdiff_cc);
[~,~,pdiff_nc2] = cyclic_map_seam(lat_g,lon_g,pdiff_nc);
[~,~,pdiff_gc2] = cyclic_map_seam(lat_g,lon_g,pdiff_gc);
[~,~,pdiff_ic2] = cyclic_map_seam(lat_g,lon_g,pdiff_ic);
[~,~,pdiff_uc2] = cyclic_map_seam(lat_g,lon_g,pdiff_uc);

[~,~,pdiff_cz2] = cyclic_map_seam(lat_g,lon_g,pdiff_cz);
[~,~,pdiff_nz2] = cyclic_map_seam(lat_g,lon_g,pdiff_nz);
[~,~,pdiff_gz2] = cyclic_map_seam(lat_g,lon_g,pdiff_gz);
[~,~,pdiff_iz2] = cyclic_map_seam(lat_g,lon_g,pdiff_iz);
[~,~,pdiff_uz2] = cyclic_map_seam(lat_g,lon_g,pdiff_uz);

%% some sort of trophic amp calc
% account for same or diff signs
czpos = find(pdiff_cz2(:)>=0);
czneg = find(pdiff_cz2(:)<0);
ccpos = find(pdiff_cc2(:)>=0);
ccneg = find(pdiff_cc2(:)<0);
cidpp = intersect(czpos,ccpos);
cidnn = intersect(czneg,ccneg);
cid1 = [cidpp;cidnn];
cid2 = intersect(czpos,ccneg);
cid3 = intersect(czneg,ccpos);

nzpos = find(pdiff_nz2(:)>=0);
nzneg = find(pdiff_nz2(:)<0);
ncpos = find(pdiff_nc2(:)>=0);
ncneg = find(pdiff_nc2(:)<0);
nidpp = intersect(nzpos,ncpos);
nidnn = intersect(nzneg,ncneg);
nid1 = [nidpp;nidnn];
nid2 = intersect(nzpos,ncneg);
nid3 = intersect(nzneg,ncpos);

gzpos = find(pdiff_gz2(:)>=0);
gzneg = find(pdiff_gz2(:)<0);
gcpos = find(pdiff_gc2(:)>=0);
gcneg = find(pdiff_gc2(:)<0);
gidpp = intersect(gzpos,gcpos);
gidnn = intersect(gzneg,gcneg);
gid1 = [gidpp;gidnn];
gid2 = intersect(gzpos,gcneg);
gid3 = intersect(gzneg,gcpos);

izpos = find(pdiff_iz2(:)>=0);
izneg = find(pdiff_iz2(:)<0);
icpos = find(pdiff_ic2(:)>=0);
icneg = find(pdiff_ic2(:)<0);
iidpp = intersect(izpos,icpos);
iidnn = intersect(izneg,icneg);
iid1 = [iidpp;iidnn];
iid2 = intersect(izpos,icneg);
iid3 = intersect(izneg,icpos);

uzpos = find(pdiff_uz2(:)>=0);
uzneg = find(pdiff_uz2(:)<0);
ucpos = find(pdiff_uc2(:)>=0);
ucneg = find(pdiff_uc2(:)<0);
uidpp = intersect(uzpos,ucpos);
uidnn = intersect(uzneg,ucneg);
uid1 = [uidpp;uidnn];
uid2 = intersect(uzpos,ucneg);
uid3 = intersect(uzneg,ucpos);

%% troph amp as ratio
cta1 = nan(size(lon_s));
cta2 = cta1;
cta3 = cta1;
nta1 = cta1;
nta2 = cta1;
nta3 = cta1;
gta1 = cta1;
gta2 = cta1;
gta3 = cta1;
ita1 = cta1;
ita2 = cta1;
ita3 = cta1;
uta1 = cta1;
uta2 = cta1;
uta3 = cta1;

%same sign
cta1(cid1) = pdiff_cz2(cid1)./pdiff_cc2(cid1);
nta1(nid1) = pdiff_nz2(nid1)./pdiff_nc2(nid1);
gta1(gid1) = pdiff_gz2(gid1)./pdiff_gc2(gid1);
ita1(iid1) = pdiff_iz2(iid1)./pdiff_ic2(iid1);
uta1(uid1) = pdiff_uz2(uid1)./pdiff_uc2(uid1);

%neg to pos
cta2(cid2) = abs(pdiff_cz2(cid2)./pdiff_cc2(cid2));
nta2(nid2) = abs(pdiff_nz2(nid2)./pdiff_nc2(nid2));
gta2(gid2) = abs(pdiff_gz2(gid2)./pdiff_gc2(gid2));
ita2(iid2) = abs(pdiff_iz2(iid2)./pdiff_ic2(iid2));
uta2(uid2) = abs(pdiff_uz2(uid2)./pdiff_uc2(uid2));

%pos to neg
cta3(cid3) = abs(pdiff_cz2(cid3)./pdiff_cc2(cid3));
nta3(nid3) = abs(pdiff_nz2(nid3)./pdiff_nc2(nid3));
gta3(gid3) = abs(pdiff_gz2(gid3)./pdiff_gc2(gid3));
ita3(iid3) = abs(pdiff_iz2(iid3)./pdiff_ic2(iid3));
uta3(uid3) = abs(pdiff_uz2(uid3)./pdiff_uc2(uid3));

%% Save

save('troph_amp_ssp585_chl_mesoz200_10yr.mat','lat_s','lon_s',...
    'cta1','nta1','gta1','ita1','uta1',...
    'cta2','nta2','gta2','ita2','uta2',...
    'cta3','nta3','gta3','ita3','uta3')

%%
load('troph_amp_ssp585_chl_mesoz200_10yr.mat')

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines;

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%%
f6 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - 1can
subplot('Position',[0.025 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta1)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Same sign','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - 1cnrm
subplot('Position',[0.025 0.6 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nta1)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - 1gfdl
subplot('Position',[0.025 0.4 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gta1)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cbr = colorbar('Position',[0.935 0.375 0.025 0.25]);
%xlabel(cbr,'ratio')

%4 - 1IPSL
subplot('Position',[0.025 0.2 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ita1)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - 1uk
subplot('Position',[0.025 0.0 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uta1)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - 2CAN
subplot('Position',[0.33 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta2)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Neg to Pos','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%7 - 2CNRM
subplot('Position',[0.33 0.6 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nta2)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - 2GFDL
subplot('Position',[0.33 0.4 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gta2)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - 2IPSL
subplot('Position',[0.33 0.2 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ita2)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%10 - 2uk
subplot('Position',[0.33 0.0 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uta2)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);


% 11 - 3CAN
subplot('Position',[0.63 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta3)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Pos to Neg','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%12 - 3CNRM
subplot('Position',[0.63 0.6 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,nta3)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'CNRM','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%13 - 3GFDL
subplot('Position',[0.63 0.4 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,gta3)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'GFDL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%14 - 3IPSL
subplot('Position',[0.63 0.2 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,ita3)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'IPSL','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%15 - 3uk
subplot('Position',[0.63 0.0 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,uta3)
colormap(cmYOR)
caxis([1 3])
text(-1.85,1.75,'UK','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_chl_zmeso_troph_amp_Hist_SSP585_90s_ratio.png'])
