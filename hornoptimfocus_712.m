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
dblim = -15;
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

dblim=-15;

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
%         i_error = loptim.field2error(apu1, ndemou1, dblim);
        apu1(apu1<dblim) = dblim;
%         i_error = immse(apu1,ndemou1)*numel(apu1);
        i_error = sum((abs(apu1)-abs(ndemou1)).^2,"all");
        all_errors(angles,i) = i_error;
        % cascade
        for j = 1:length(coeffs)
            new_coeffs = coeffs;
            new_coeffs(j) = new_coeffs(j) + lr*randi([-1,1])*rand(1,1);
%             new_coeffs(j) = new_coeffs(j) + lr*rand(1,1);
            aperthickness =l.phaseprofile(new_coeffs, antenna_r);
            aperlens = l.makelensfromthickness(aperthickness, antenna_r, 1);
            apu1 = l.lenspropagate(u0, aperlens, 0, zmid);
            apu1 = l.lenspropagate(apu1, focuslens, 0, z2);
%             j_error = loptim.field2error(apu1, ndemou1, dblim);
            apu1(apu1<dblim) = dblim;
%             j_error =  immse(apu1,ndemou1)*numel(apu1);
            j_error = sum((abs(apu1)-abs(ndemou1)).^2,"all");
            % currently just a one step thing
            % can implement simulated annealing later
%             disp(del_error);
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

%% Overall error

% Plot accumulated error
figure('Name', 'Plot of Error')
plot(sum(all_errors),"-","LineWidth",3);
title("Error for Optimizing over Iterations")
xlabel("Iteration");
ylabel("MSE");
grid on;
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1, 'FontWeight','Normal');


%% Error Plot

num=5;

% Plot accumulated error
figure('Name', 'Plot of Error')
plot(all_errors(num,:),"-","LineWidth",3,"Color",'#7E2F8E');
% plot(sum(all_errors)/length(xangles),"-x","LineWidth",1);
title(strjoin(["Error for optimizing", int2str(xangles(num)),"degrees"]));
% title("Error for Optimizing over Iterations")
xlabel("Iteration");
ylabel("MSE");
grid on;
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1, 'FontWeight','Normal');

%% playground

a = [1:10 ; 11:20; 21:30; 31:40; 41:50];
b = reshape(a.', 1,50);

%% Error together without summing 

% Plot accumulated error
figure('Name', 'Plot of Error');
all_errors_all = reshape(all_errors.',1,500);
plot(all_errors_all,"LineWidth",3,"Color",'#7E2F8E');
hold on;
for i = 1:4
    xline(i*100,"--","LineWidth",2);
end
hold off;
% plot(sum(all_errors)/length(xangles),"-x","LineWidth",1);
title(strjoin(["Error for optimizing over all iterations"]));
% title("Error for Optimizing over Iterations")
xlabel("Iteration");
ylabel("MSE");
grid on;
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1, 'FontWeight','Normal');

%% Accumulated Error Plot
errors_over_sep = [
2.142894554261546,
2.016048991405506,
1.804485864787087,
2.263408667827682,
1.339249745141558,
2.008709000466991,
1.775916624046341,
2.461041981788260,
2.303938072130269,
2.379091618931086,
];

separation_array = 1:10;


figure("Name","Weighted Error over Separation")
plot(separation_array, errors_over_sep,'--x','LineWidth',5,'MarkerSize',10,'Color','#7E2F8E','LineWidth',1.5);
grid on;
title("Error for different separation values");
yline(min(errors_over_sep),'--');
xlim([0 11]);
xlabel("Separation distance (mm)");
ylabel("Error");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1, 'FontWeight','Normal');



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

dblim = -25;

hornangle = 50;

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

% figure;
% subplot(121);
% imagesc(p.x, p.y, ncpu1);
% title(["Regular Lens at ", int2str(hornangle)," degrees"]);
% caxis([-20 0]);
% axis square;
% colormap jet;
% colorbar;
% 
% subplot(122);
% imagesc(p.x, p.y, napu1);
% title(["Aperture Lens at ", int2str(hornangle)," degrees"]);
% axis square;
% colormap jet;
% caxis([-20 0]);
% colorbar;

figure;
imagesc(p.x, p.y, napu1);
title(["Aperture Lens at ", int2str(hornangle)," degrees"]);
axis square;
colormap jet;
caxis([-20 0]);
colorbar;
set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1, 'FontWeight','Normal');



figure('Name', 'Cross section');
plot(p.x, ncpu1(ceil(M/2),:));
hold on;
plot(p.x, napu1csection);
% plot(p.x, ndemou1(ceil(M/2),:));
hold off;
legend("cp lens", "aper lens");
title("Cross section intensity profile");

%% Ideal field distribution

figure;
imagesc(p.x, p.y, ndemou1);
caxis([dblim 0]);
axis square;
colormap jet;
colorbar;

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
