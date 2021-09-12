close all;
format long;
addpath('./Classes');


lambda = 1e-3;
k0 = 2*pi/lambda;
antenna_r = 10e-3;
dx = 0.5*lambda;
L = 200e-3;
M = L/dx + 1;
z1 = 56e-3;
z2 = 40e-3;

k = -824;
R = 1/0.128;

% Instantiate classes
p = Propagator(L, lambda, dx);
l = Lens(L, lambda, dx);

% quick lambda function
normdb = @(u1) mag2db(abs(u1)) - max(max(mag2db(abs(u1))));

coeffs = [-13.5078 22.5578 -5.3161 6.8526 -2.1028];
r1 = p.x/antenna_r;
aperlens = l.makephaselens(coeffs, antenna_r);
cplens = l.makecplens(z1, z2, antenna_r);

u0 = p.pso(0, 0, z1);
u1 = l.lenspropagate(u0, aperlens, z1, z2);

figure;
imagesc(p.x, p.y, normdb(u1));
caxis([-15 0]);
colorbar;
axis square;
colormap jet;


    