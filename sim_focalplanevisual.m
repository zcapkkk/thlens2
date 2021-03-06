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
zmid = 5e-3;
dblim = -20;

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



% initialize coefficients

ape = readmatrix('Lenses/1_aperturelens.xlsx');
foc = readmatrix('Lenses/1_focuslens.xlsx');

% prepare the lens plane
% pd0 = (length(xx)-length(ape))/2;
% aperlens = padarray(ape, [pd0 pd0], 1, 'both');
% focuslens = padarray(foc,[pd0 pd0], 1, 'both');
coeffs = [-13.257440755047092  22.817438177933152...
    -5.072199588228340 7.110879608391324  -1.859533046641054];
aperlens = l.makephaselens(coeffs, antenna_r, 1);
focuslens = l.makecplens(z1+zmid, z2, antenna_r, 1);

%% Propagation

for angles = 1:length(xangles)    
    pointsource = p.pso(xangles(angles), 0, z1);
    u0_phase = angle(p.prop(pointsource, z1*cosd(xangles(angles))));

    ampdata = readmatrix(int2str(xangles(angles))+"deg_Amp");
    phasedata = readmatrix(int2str(xangles(angles))+"deg_Phase");
    phasedata = deg2rad(phasedata);
    u0_amp = interp2(Hornx, Horny, ampdata, xx, yy, 'nearest', 0);
    u0 = u0_amp.*exp(1i.*u0_phase);
    
    % Lenses
    
    % regular lens
    cplens = l.makecplens(z1, z2, antenna_r, 1);
    cpu1 = l.lenspropagate(u0, cplens, 0, z2);
    ncpu1 = normdb(cpu1);
    ncpu1csection = ncpu1(ceil(M/2),:);
    
    
    
    % aperture lens
    apu1 = l.lenspropagate(u0, aperlens, 0, z2);
    napu1 = normdb(apu1);
    napu1csection = napu1(ceil(M/2),:);
    
    % doublet lens
    dapu1 = l.lenspropagate(u0, aperlens, 0, zmid);
    dapu1 = l.lenspropagate(dapu1, focuslens, 0, z2);
    ndapu1 = normdb(dapu1);
    ndapu1csection = ndapu1(ceil(M/2),:);
    
    name4fig = int2str(xangles(angles)) + " degrees";
    figure("Name", name4fig);

    
    t = tiledlayout(1,2, 'TileSpacing', 'tight', 'Padding', 'tight');   
    nexttile;
    imagesc(p.x/1e-3,p.y/1e-3,imrotate(ncpu1,90));
    [t,s] = title("Conventional", ["Parabolic, ", int2str(xangles(angles))+ " degrees"]);
    colormap jet;
    c1 = colorbar;
    c1.Label.String = "dB";
    axis equal;
    ylim([-50 10]);
    xlim([-10 10]);
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
    caxis([dblim 0]);
    
    nexttile;
    imagesc(p.x/1e-3,p.y/1e-3,imrotate(ndapu1,90));
    [t,s] = title("Doublet", ["Polynomial + Parabolic", int2str(xangles(angles))+ " degrees"]);
    colormap jet;
    c3 = colorbar;
    c3.Label.String = "dB";
    axis equal;
    axis xy;
    ylim([-50 10]);
    xlim([-10 10]);
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
    caxis([dblim 0]);
 
end