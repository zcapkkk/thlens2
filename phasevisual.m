% use for visualizing phase distribution o reflectarray
close all;
addpath('./Lenses');
addpath('./Classes');
addpath('./MeasData');


L = 200e-3;
lambda = 1e-3;
dx = (1/3)*lambda;

p = Propagator(L, lambda, dx);
l = Lens(L, lambda, dx);



ape = readmatrix('Lenses/1_aperturelens.xlsx');
foc = readmatrix('Lenses/1_focuslens.xlsx');
x_a = p.x(abs(p.x)<=0.011);
y_a = p.y(abs(p.y)<=0.011);
figure;
imagesc(x_a/lambda, y_a/lambda, ape);
xlabel('x (mm)');
ylabel('y (mm)');
axis square;
colormap jet;
% title("Aperture Transmitarray", "FontWeight","Normal");
cape = colorbar;
cape.Label.String = 'Phase shift (degrees)';
set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1, 'FontWeight','Normal');
saveas(gcf, 'MeasData/aperture_phase.png');

figure;
imagesc(x_a/lambda, y_a/lambda, foc);
xlabel('x (mm)');
ylabel('y (mm)');
axis square;
colormap jet;
% title("Focusing Transmitarray", "FontWeight","Normal");
cfoc = colorbar;
cfoc.Label.String = 'Phase shift (degrees)';
set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1);
saveas(gcf, 'MeasData/focusing_phase.png');

%% Make regular lens
lambda = 1e-3;
k0 = 2*pi/lambda;
antenna_r = 10e-3;
dx = (1/3)*lambda;
L = 200e-3;
M = round(L/dx) + 1;
z1 = 56e-3;
z2 = 40e-3;
cplens = l.makecplens(z1, z2, antenna_r, 1);
cplens = angle(cplens);
cplens(cplens < 0) = 2*pi + cplens(cplens<0);
cplens = rad2deg(cplens);
figure;
imagesc(p.x/lambda, p.y/lambda, cplens);
xlabel('x (mm)');
ylabel('y (mm)');
xlim([-10 10]);
ylim([-10 10]);
axis square;
colormap jet;
% title("Conventional Transmitarray", "FontWeight","Normal");
cfoc = colorbar;
cfoc.Label.String = 'Phase shift (degrees)';
set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1);
saveas(gcf, 'MeasData/conventional_phase.png');
