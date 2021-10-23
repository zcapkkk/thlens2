% use for visualizing phase distribution o reflectarray
close all;
addpath('./Lenses');
addpath('./Classes');


L = 200e-3;
lambda = 1e-3;
dx = (1/3)*lambda;

p = Propagator(L, lambda, dx);
l = Lens(L, lambda, dx);



ape = readmatrix('Lenses/1_aperturelens.xlsx');
foc = readmatrix('Lenses/1_focuslens.xlsx');

figure;
imagesc(p.x/lambda, p.y/lambda, ape);
xlabel('mm');
ylabel('mm');
axis square;
colormap jet;
title("Aperture Transmitarray");
cape = colorbar;
cape.Label.String = 'Phase shift (degrees)';

figure;
imagesc(p.x/lambda, p.y/lambda, foc);
xlabel('mm');
ylabel('mm');
axis square;
colormap jet;
title("Focusing Transmitarray");
cfoc = colorbar;
cfoc.Label.String = 'Phase shift (degrees)';