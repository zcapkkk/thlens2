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


xangles = [0 30:10:50];



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

% arrays to store 1D cross sectional results
regular = zeros(length(xangles), L/dx+1);
doublet = zeros(length(xangles), L/dx+1);

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
    regular(angles, :) = ncpu1csection;
        
    % doublet lens
    dapu1 = l.lenspropagate(u0, aperlens, 0, zmid);
    dapu1 = l.lenspropagate(dapu1, focuslens, 0, z2);
    ndapu1 = normdb(dapu1);
    ndapu1csection = ndapu1(ceil(M/2),:);
    doublet(angles, :) = ndapu1csection;
end

figure("Name","Regular");
hold on;
for i = 1:length(xangles)
    plot(p.x/lambda,regular(i,:),LineWidth=2);
end
% xline([(z2)*tand(xangles)/lambda],'--','LineWidth',1.5);
hold off;
grid MINOR;
xlim([-6 60]);
ylim([-20 0]);
legend(["$0^\circ$","$30^\circ$","$40^\circ$",...
    "$50^\circ$"],"Interpreter","latex");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel("y (mm)");
ylabel("dB");

figure("Name","Doublet");
hold on;
for i = 1:length(xangles)
    plot(p.x/lambda,doublet(i,:),'LineWidth',2);
end
% xline([(z2)*tand(xangles)/lambda],'--','LineWidth',1.5);
hold off;
grid MINOR;
xlim([-6 60]);
ylim([-20 0]);
legend(["$0^\circ$","$30^\circ$","$40^\circ$",...
    "$50^\circ$"],"Interpreter","latex");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel("y (mm)");
ylabel("dB");