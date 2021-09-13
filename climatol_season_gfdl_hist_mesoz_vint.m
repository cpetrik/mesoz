% Calc seasonal climatology
% Test by hand vs. ClimateDataToolbox

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';

load([fpath 'gfdl_hist_zmeso_vint_monthly_1950_2014.mat'])

%% ClimateDataToolbox
runs = find(yr>1950 & yr<=2015);
%zmc = climatology(zmeso_vint,runs);
[zmc,tc] = climatology(zmeso_vint,runs,'monthly');

%% By hand
nmo = length(runs);
nyr = length(runs)/12;

zm_mo = nan*ones(360,180,12);
for m = 1:12
    mo = m:12:nyr;
    zm_mo(:,:,m) = nanmean(double(zmeso_vint(:,:,mo)),3);
end

%% Climatology package
figure

h = imagescn(lon',lat',real(log10(zmc(:,:,1)))');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title(datestr(datenum(0,0,1),'mmmm'))
caxis([-1.5 0])

hold on
he = earthimage;
uistack(he,'bottom')

gif('gfdl_mesoz_monthly_climatology.gif','frame',gcf,'delaytime',1/12,'nodither')

for k=2:12
    h.CData = real(log10(zmc(:,:,k)))';
    title(datestr(datenum(0,k,1),'mmmm'))
    gif
end

%% By hand
figure

h = imagescn(lon',lat',real(log10(zm_mo(:,:,1)))');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title(datestr(datenum(0,0,1),'mmmm'))
caxis([-1.5 0])

hold on
he = earthimage;
uistack(he,'bottom')

gif('gfdl_mesoz_monthly_climatol_calc.gif','frame',gcf,'delaytime',1/12,'nodither')

for k=2:12
    h.CData = real(log10(zm_mo(:,:,k)))';
    title(datestr(datenum(0,k,1),'mmmm'))
    gif
end

%% By hand shows seasonal progression between hemispheres
zmc_DJF = mean(zm_mo(:,:,[1 2 12]),3);
zmc_MAM = mean(zm_mo(:,:,3:5),3);
zmc_JJA = mean(zm_mo(:,:,6:8),3);
zmc_SON = mean(zm_mo(:,:,9:11),3);
zmc_all = mean(zm_mo,3);

%%
figure(2)
h = imagescn(lon',lat',log10(zmc_DJF)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title('Winter')
caxis([-1.5 0])
hold on
he = earthimage;
uistack(he,'bottom')

figure(3)
h = imagescn(lon',lat',log10(zmc_MAM)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title('Spring')
caxis([-1.5 0])
hold on
he = earthimage;
uistack(he,'bottom')

figure(4)
h = imagescn(lon',lat',log10(zmc_JJA)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title('Summer')
caxis([-1.5 0])
hold on
he = earthimage;
uistack(he,'bottom')

figure(5)
h = imagescn(lon',lat',log10(zmc_SON)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean algae
title('Fall')
caxis([-1.5 0])
hold on
he = earthimage;
uistack(he,'bottom')

%%
figure
h = imagescn(lon',lat',(zmc_DJF - zmc_JJA)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean balance
title('Winter - Summer')
caxis([-0.03 0.03])
hold on
he = earthimage;
uistack(he,'bottom')

figure
h = imagescn(lon',lat',(zmc_DJF - zmc_MAM)');
cb = colorbar;
ylabel(cb,['mesoz abund (' units ')'])
cmocean balance
title('Winter - Spring')
caxis([-0.03 0.03])
hold on
he = earthimage;
uistack(he,'bottom')

%%
save([fpath 'gfdl_hist_zmeso_vint_climatol_1950_2014.mat'],'yr',...
    'long_name','standard_name','units','lat','lon','runs',...
    'zmc_DJF','zmc_MAM','zmc_JJA','zmc_SON','zmc_all');

