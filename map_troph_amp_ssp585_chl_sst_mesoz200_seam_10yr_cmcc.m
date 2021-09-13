% CMIP6 output 
% mesoz 200m integrations
% Map ESM future diffs 2090s vs. 1990s

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/Fish-MIP/CMIP6/driver_analysis/zmeso_figs/';

load('cmip6_hist_ssp585_space_means_90s_zmeso200_schl_sst_same_orientation.mat')

%% Diffs 
diff_ct = FCsst90 - HCsst90;
diff_mt = FMsst90 - HMsst90;
diff_nt = FNsst90 - HNsst90;
diff_gt = FGsst90 - HGsst90;
diff_it = FIsst90 - HIsst90;
diff_ut = FUsst90 - HUsst90;

%% Pdiffs
pdiff_cc = (FCchl - HCchl) ./ HCchl;
pdiff_mc = (FMchl - HMchl) ./ HMchl;
pdiff_nc = (FNchl - HNchl) ./ HNchl;
pdiff_gc = (FGchl - HGchl) ./ HGchl;
pdiff_ic = (FIchl - HIchl) ./ HIchl;
pdiff_uc = (FUchl - HUchl) ./ HUchl;

pdiff_cz = (FCzm - HCzm) ./ HCzm;
pdiff_mz = (FMzm - HMzm) ./ HMzm;
pdiff_nz = (FNzm - HNzm) ./ HNzm;
pdiff_gz = (FGzm - HGzm) ./ HGzm;
pdiff_iz = (FIzm - HIzm) ./ HIzm;
pdiff_uz = (FUzm - HUzm) ./ HUzm;

%% Fix seam
[lat_s,lon_s,diff_ct2] = cyclic_map_seam(lat_g,lon_g,diff_ct);
[~,~,diff_mt2] = cyclic_map_seam(lat_g,lon_g,diff_mt);
[~,~,diff_nt2] = cyclic_map_seam(lat_g,lon_g,diff_nt);
[~,~,diff_gt2] = cyclic_map_seam(lat_g,lon_g,diff_gt);
[~,~,diff_it2] = cyclic_map_seam(lat_g,lon_g,diff_it);
[~,~,diff_ut2] = cyclic_map_seam(lat_g,lon_g,diff_ut);

[~,~,pdiff_cc2] = cyclic_map_seam(lat_g,lon_g,pdiff_cc);
[~,~,pdiff_mc2] = cyclic_map_seam(lat_g,lon_g,pdiff_mc);
[~,~,pdiff_nc2] = cyclic_map_seam(lat_g,lon_g,pdiff_nc);
[~,~,pdiff_gc2] = cyclic_map_seam(lat_g,lon_g,pdiff_gc);
[~,~,pdiff_ic2] = cyclic_map_seam(lat_g,lon_g,pdiff_ic);
[~,~,pdiff_uc2] = cyclic_map_seam(lat_g,lon_g,pdiff_uc);

[~,~,pdiff_cz2] = cyclic_map_seam(lat_g,lon_g,pdiff_cz);
[~,~,pdiff_mz2] = cyclic_map_seam(lat_g,lon_g,pdiff_mz);
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

mzpos = find(pdiff_mz2(:)>=0);
mzneg = find(pdiff_mz2(:)<0);
mcpos = find(pdiff_mc2(:)>=0);
mcneg = find(pdiff_mc2(:)<0);
midpp = intersect(mzpos,mcpos);
midnn = intersect(mzneg,mcneg);
mid1 = [midpp;midnn];
mid2 = intersect(mzpos,mcneg);
mid3 = intersect(mzneg,mcpos);

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
mta1 = cta1;
mta2 = cta1;
mta3 = cta1;
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
mta1(mid1) = pdiff_mz2(mid1)./pdiff_mc2(mid1);
nta1(nid1) = pdiff_nz2(nid1)./pdiff_nc2(nid1);
gta1(gid1) = pdiff_gz2(gid1)./pdiff_gc2(gid1);
ita1(iid1) = pdiff_iz2(iid1)./pdiff_ic2(iid1);
uta1(uid1) = pdiff_uz2(uid1)./pdiff_uc2(uid1);

%neg to pos
cta2(cid2) = abs(pdiff_cz2(cid2)./pdiff_cc2(cid2));
mta2(mid2) = abs(pdiff_mz2(mid2)./pdiff_mc2(mid2));
nta2(nid2) = abs(pdiff_nz2(nid2)./pdiff_nc2(nid2));
gta2(gid2) = abs(pdiff_gz2(gid2)./pdiff_gc2(gid2));
ita2(iid2) = abs(pdiff_iz2(iid2)./pdiff_ic2(iid2));
uta2(uid2) = abs(pdiff_uz2(uid2)./pdiff_uc2(uid2));

%pos to neg
cta3(cid3) = abs(pdiff_cz2(cid3)./pdiff_cc2(cid3));
mta3(mid3) = abs(pdiff_mz2(mid3)./pdiff_mc2(mid3));
nta3(nid3) = abs(pdiff_nz2(nid3)./pdiff_nc2(nid3));
gta3(gid3) = abs(pdiff_gz2(gid3)./pdiff_gc2(gid3));
ita3(iid3) = abs(pdiff_iz2(iid3)./pdiff_ic2(iid3));
uta3(uid3) = abs(pdiff_uz2(uid3)./pdiff_uc2(uid3));

%% Save

save('troph_amp_ssp585_chl_mesoz200_10yr.mat','lat_s','lon_s',...
    'cta1','mta1','nta1','gta1','ita1','uta1',...
    'cta2','mta2','nta2','gta2','ita2','uta2',...
    'cta3','mta3','nta3','gta3','ita3','uta3')

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

%1 - 1cmcc
subplot('Position',[0.025 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mta1)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Same sign','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CMCC','HorizontalAlignment','center')
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

% 6 - 2Cmcc
subplot('Position',[0.33 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mta2)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Neg to Pos','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CMCC','HorizontalAlignment','center')
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


% 11 - 3Cmcc
subplot('Position',[0.63 0.8 0.3 0.2])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,mta3)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Pos to Neg','HorizontalAlignment','center','FontWeight','bold')
text(-1.85,1.75,'CMCC','HorizontalAlignment','center')
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

%% CAN only
figure
%1 - 1can
subplot('Position',[0.015 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta1)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Same sign','HorizontalAlignment','center','FontWeight','bold')
%text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 6 - 2CAN
subplot('Position',[0.32 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta2)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Neg to Pos','HorizontalAlignment','center','FontWeight','bold')
%text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 11 - 3CAN
subplot('Position',[0.625 0.5 0.3 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat_s,lon_s,cta3)
colormap(cmYOR)
caxis([1 3])
text(0,1.75,'Pos to Neg','HorizontalAlignment','center','FontWeight','bold')
%text(-1.85,1.75,'CAN','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cbr = colorbar('Position',[0.935 0.5 0.025 0.3]);

print('-dpng',[ppath 'Map_chl_zmeso_troph_amp_Hist_SSP585_90s_ratio_CANonly.png'])
