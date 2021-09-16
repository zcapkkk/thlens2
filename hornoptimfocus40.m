close all;
format long;
addpath('./Classes');
addpath('./HornData/Dist56');


lambda = 1e-3;
k0 = 2*pi/lambda;
antenna_r = 10e-3;
dx = (1/3)*lambda;
L = 200e-3;
M = round(L/dx) + 1;
z1 = 56e-3;
z2 = 40e-3;
dblim = -25;
zmid = 5e-3;

% Instantiate classes
p = Propagator(L, lambda, dx);
l = Lens(L, lambda, dx);
loptim = LensOptimizer(L, lambda, dx, z1, z2, antenna_r);

% quick lambda function
normdb = @(u1) mag2db(abs(u1)) - max(max(mag2db(abs(u1))));


xangles = [0 30:10:60];



% horn setup
hornx = 1e-3*(-10:0.2:10);
horny = 1e-3*(-10:0.2:10);
[Hornx, Horny] = meshgrid(hornx, horny);
[xx, yy] = meshgrid(p.x, p.y);


% optim setup:



% % initialize coefficients
% coeffs =  [-13.257589733620399, 22.810501285918946, -5.075819033570514,...
%     7.082499034482761, -1.882158081684554];
coeffs = [-13.257440755047092  22.817438177933152  -5.072199588228340 7.110879608391324  -1.859533046641054];
aperlens = l.makephaselens(coeffs, antenna_r, 1); 
focusthickness = l.cpthickness(z1+zmid, z2);
focuslens = l.makelensfromthickness(focusthickness, antenna_r, 1);

%% Testing the lens

hornangle = 40;

% make source
pointsource = p.pso(hornangle, 0, z1);
u0_phase = angle(p.prop(pointsource, z1*cosd(hornangle)));


ampdata = readmatrix(int2str(hornangle)+"deg_Amp");
phasedata = readmatrix(int2str(hornangle)+"deg_Phase");
phasedata = deg2rad(phasedata);
u0_amp = interp2(Hornx, Horny, ampdata, xx, yy, 'nearest', 0);
u0 = u0_amp.*exp(1i.*u0_phase);

% aperture lens
apu1 = l.lenspropagate(u0, aperlens, 0, zmid);
apu1 = l.lenspropagate(apu1, focuslens, 0, z2);
napu1 = normdb(apu1);
napu1csection = napu1(ceil(M/2),:);

% Get aperture lens data
bapu1 = imbinarize(napu1, -5);
binfo = regionprops(bapu1, "Eccentricity","MajorAxisLength","MinorAxisLength");
disp(binfo.Eccentricity)
disp(binfo.MajorAxisLength)
disp(binfo.MinorAxisLength)


figure;
imagesc(p.x, p.y, bapu1);
title(["Aperture Lens at ", int2str(hornangle)," degrees"]);
axis square;
colormap jet;
colorbar;


