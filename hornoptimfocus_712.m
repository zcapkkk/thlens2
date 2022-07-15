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


xangles = [0 30:10:69];



% horn setup
hornx = 1e-3*(-10:0.2:10);
horny = 1e-3*(-10:0.2:10);
[Hornx, Horny] = meshgrid(hornx, horny);
[xx, yy] = meshgrid(p.x, p.y);


% coeffs =  [-13.257589733620399, 22.810501285918946, -5.075819033570514,...
%     7.082499034482761, -1.882158081684554];
focuslens = l.makecplens(z1+zmid, z2, antenna_r, 1);

%% training


steps = 100;
lr = 1;
prev_error = 1e6;
coeffs = [0 0 0 0 0];
% coeffs = randi(5,1,5)-randi(10,1,5);
disp(coeffs);
all_errors = zeros(length(xangles),steps);
tic 
disp("Begin optimization...")
for angles = 1:length(xangles)
    disp(xangles(angles));
    % ideal focus
    demou1 = loptim.idealcpu1(xangles(angles),0);
    ndemou1 = normdb(demou1);
    ndemou1csection = ndemou1(ceil(M/2),:);
    ndemou1(ndemou1 < dblim) = dblim;
    pointsource = p.pso(xangles(angles), 0, z1);
    u0_phase = angle(p.prop(pointsource, z1*cosd(xangles(angles))));


    ampdata = readmatrix(int2str(xangles(angles))+"deg_Amp");
%     phasedata = readmatrix(int2str(xangles(angles))+"deg_Phase");
%     phasedata = deg2rad(phasedata);
    u0_amp = interp2(Hornx, Horny, ampdata, xx, yy, 'nearest', 0);
    u0 = u0_amp.*exp(1i.*u0_phase);

    
    for i = 1:steps
        T = steps - i;
        aperthickness = l.phaseprofile(coeffs, antenna_r);
        aperlens = l.makelensfromthickness(aperthickness, antenna_r, 1);
        apu1 = l.lenspropagate(u0, aperlens, 0, zmid);
        apu1 = l.lenspropagate(apu1, focuslens, 0, z2);
        i_error = loptim.field2error(apu1, ndemou1, dblim);
        all_errors(angles,i) = i_error;
        % cascade
        for j = 1:length(coeffs)
            new_coeffs = coeffs;
            new_coeffs(j) = new_coeffs(j) + lr*randi([-1,1])*rand(1,1);
            aperthickness =l.phaseprofile(new_coeffs, antenna_r);
            aperlens = l.makelensfromthickness(aperthickness, antenna_r, 1);
            apu1 = l.lenspropagate(u0, aperlens, 0, zmid);
            apu1 = l.lenspropagate(apu1, focuslens, 0, z2);
            j_error = loptim.field2error(apu1, ndemou1, dblim);
            % currently just a one step thing
            % can implement simulated annealing later
            del_error = j_error - i_error;
            if del_error < 0
                coeffs = new_coeffs;
            else
                prob = exp(-del_error/T);
                if rand(1,1) < prob
                    coeffs = new_coeffs;
                end
            end
        end
    end
%     disp([i, coeffs])          
end

toc

aperthickness =l.phaseprofile(coeffs, antenna_r);
aperlens = l.makelensfromthickness(aperthickness, antenna_r, 1);

disp(coeffs)
%% Error Plot

num=1;

% Plot accumulated error
figure('Name', 'Plot of Error')
% plot(all_errors(num,:),"-x","LineWidth",1,"Color",'#7E2F8E');
plot(sum(all_errors)/length(xangles),"-x","LineWidth",1);
% title(strjoin(["Error for optimizing", int2str(xangles(num)),"degrees"]));
title("Error for Optimizing over Iterations")
xlabel("Iteration");
ylabel("MSE");
grid on;
set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1, 'FontWeight','Normal');
%% Plotting

figure('Name', 'Phase profile');
subplot(131);
pf = 0;
for i = 1:length(coeffs)
    pf = pf + coeffs(i)*(p.x/antenna_r).^(2*i);
end
plot(p.x, pf);
xlim([-antenna_r antenna_r]);
axis square;
subplot(132);
focuslensphase = angle(aperlens);
plot(p.x, focuslensphase(ceil(M/2),:));
xlim([-antenna_r antenna_r]);
axis square;
subplot(133);
imagesc(p.x, p.y, focuslensphase);
xlim([-antenna_r antenna_r]);
ylim([-antenna_r antenna_r]);
title("Lens phase dist.");
axis square;
colormap jet;
colorbar;


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

% propagate through aperture lens
apu1 = l.lenspropagate(u0, aperlens, 0, zmid);
apu1 = l.lenspropagate(apu1, focuslens, 0, z2);
napu1 = normdb(apu1);
napu1(napu1 < dblim) = dblim;
napu1csection = napu1(ceil(M/2),:);


% regular lens
cplens = l.makecplens(z1, z2+zmid, antenna_r, 1);
cpu1 = l.lenspropagate(u0, cplens, 0, z2);
ncpu1 = normdb(cpu1);
ncpu1csection = ncpu1(ceil(M/2),:);
ncpu1(ncpu1 < dblim) = dblim;
% benchmarks
cpdbtotal = sum(ncpu1, 'all');
cperror = immse(ndemou1, ncpu1);

% figure;
% imagesc(p.x, p.y, normdb(demou1));
% caxis([-20 0]);
% axis square;
% colormap jet;
% colorbar;


% do NOT use normalised to show absolute loss

figure;
subplot(121);
imagesc(p.x, p.y, ncpu1);
title(["Regular Lens at ", int2str(hornangle)," degrees"]);
caxis([-20 0]);
axis square;
colormap jet;
colorbar;

subplot(122);
imagesc(p.x, p.y, napu1);
title(["Aperture Lens at ", int2str(hornangle)," degrees"]);
axis square;
colormap jet;
caxis([-20 0]);
colorbar;

figure('Name', 'Cross section');
plot(p.x, ncpu1(ceil(M/2),:));
hold on;
plot(p.x, napu1csection);
% plot(p.x, ndemou1(ceil(M/2),:));
hold off;
legend("cp lens", "aper lens");
title("Cross section intensity profile");

%% Export the lens
% disp(coeffs)
% 
% st = ceil(M/2) - round(antenna_r/dx);
% ed = ceil(M/2) + round(antenna_r/dx);
%  
% buildfocuslensphase = rad2deg(focuslensphase);
% buildfocuslensphase = buildfocuslensphase + 180;
% buildfocuslensphase(p.X.^2+p.Y.^2>antenna_r^2) = 0;
% buildfocuslensphase = buildfocuslensphase(st:ed, st:ed);
% delete 'Lenses/focuslens2.xlsx';
% writematrix(buildfocuslensphase, 'Lenses/focuslens2.xlsx');
