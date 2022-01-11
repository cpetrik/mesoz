% Orders of mag diff in obs
% Collated in Stock et al. 2014 Tables 7-9

clear all
close all

zmeso(1) = 286;     %HOT
zmeso(2) = 137;     %BATS
zmeso(3) = 1320;    %JGOFS Arabian Sea
zmeso(4) = 350;     %JGOFS EEQPAC
zmeso(5) = 2125;    %PAPA
zmeso(6) = 88;      %NAtl Bloom 1
zmeso(7) = 2190;    %NAtl Bloom 2

%%
zmax = max(zmeso);
zmin = min(zmeso);
log10(zmax) - log10(zmin);

%% COPEPOD
load('copepod-2012_cmass_all_gridded.mat','lat_g','lon_g',...
    'zoo_g','fileid','units')

% From m-3 to m-2 (integrate top 200m)
zoo_200 = zoo_g*200;

zmax = nanmax(zoo_g(:));
zmin = nanmin(zoo_g(:));
log10(zmax) - log10(zmin)