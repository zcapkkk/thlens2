%Calculate the radiation pattern 
clear all;close all;
clc;

deg00 = "20211230_B00deg.mat";
deg20 = "20220102_B20deg.mat";
deg30 = "20211231_B30deg.mat";
deg40 = "20211230_B40deg.mat";
deg45 = "20220102_B45deg.mat";
deg50 = "20211231_B50deg.mat";

degrees = [deg00, deg20, deg30, deg40, deg45, deg50];
big_mat = zeros(9, length(degrees));

for ii = 1:length(degrees)
    fname = convertStringsToChars(degrees(ii));
    load(fname);
    for NumPoint = 1:9
        sz=size(sdata.s21);
        M=sz(1);
        N=sz(2);
        MagEx=ones(M,N);
        PhaseEx=ones(M,N);
        for m=1:M
            for n=1:N
                Temp=cell2mat(sdata.s21(m,n));
                MagEx(m,n)=abs(Temp(NumPoint,1));
                PhaseEx(m,n)=angle(Temp(NumPoint,1))/pi*180;
            end
        end
        Ex_inc=MagEx.*exp(1i*PhaseEx/180*pi);
        MagEx_dB=20*log10(MagEx);
        MagEx=MagEx/max(MagEx(:));
        big_mat(NumPoint,ii) = max(max(MagEx_dB));
    end
end

 
%% Scanning Loss Plotting

% big_mat = [-26.1956 ,   -22.7192 ,   -25.2628 ,   -27.7436 ,   -28.2097 ,   -31.4879; 
%   -25.0033 ,   -24.9864 ,   -26.3667 ,   -27.6185 ,   -29.7617 ,   -32.4617;
% -26.4165 ,   -27.1177 ,   -27.4988 ,   -30.0747 ,   -32.3394 ,   -35.2886;
%  -21.6082 ,   -18.7023 ,   -22.4338 ,   -24.1148 ,   -25.0187 ,   -28.3977;
%  -21.6666 ,   -22.4729 ,   -23.1229 ,   -25.9072 ,   -27.8629 ,   -30.9582;
% -23.9707 ,   -20.4371 ,   -22.1450 ,   -24.6549 ,   -27.2830 ,   -30.0080;
%  -23.9625 ,   -22.7944 ,   -25.3644 ,   -27.7386 ,   -28.8239 ,   -31.8823];
big_mat = big_mat - big_mat(:,1);
x_degrees=[0,20:10:40,45,50];
tscanloss = [8.366140278832628   6.956136084256621   6.210060489485697   3.769371701218659   0.082259739108424];
tscanloss = tscanloss(1:end-1);
legend_big_mat = 240:10:320;
% zero_ave = mean(big_mat,1);
% zero_ave = zero_ave(1);
ghz_strjoin = @(x) strjoin([int2str(x),"GHz"],'');
legend_big_mat = arrayfun(ghz_strjoin,legend_big_mat);
figure;
hold on;
for i = 1:length(big_mat)
    plot(x_degrees,big_mat(i,:),'--x',LineWidth=2,MarkerSize=15);
end
plot([0,30:10:50],tscanloss-max(tscanloss),'-o',LineWidth=5,MarkerSize=15,Color=[0 0.4470 0.7410 0.5]);
legend([legend_big_mat,"Theory"],Location="southwest",NumColumns=2);
caxis([-36 -15]);
grid on;
title("Loss");
xlabel("Angle (degrees)");
ylabel("Loss (dB)");
hold off;
alpha(.5);
set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1);

    
%% Scanning Loss over frequency

sl_bm = big_mat(:,1) - big_mat(:,end);

figure;
plot(240:10:320, sl_bm,'-x',LineWidth=3,MarkerSize=15,Color='#D95319');
grid on;
title("Scanning loss (angle max to min)");
xlabel("Frequency");
ylabel("Scanning loss (dB)");
hold off;
set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1);

sl_bm2 = max(big_mat,[],2)-min(big_mat,[],2);
figure;
plot(240:10:320, sl_bm2,'-x',LineWidth=3,MarkerSize=15,Color='#0072BD');
grid on;
title("Scanning loss (loss max to min)");
xlabel("Frequency");
ylabel("Scanning loss (dB)");
hold off;
set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1);

