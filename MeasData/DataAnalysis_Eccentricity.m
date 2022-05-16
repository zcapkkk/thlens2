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
eccen = zeros(8,6);
maj_ax = zeros(8,6);
min_ax = zeros(8,6);


for ii = 1:length(degrees)
    fname = convertStringsToChars(degrees(ii));
    load(fname);
    for jj = 1:8
    

        sz=size(sdata.s21);
        M=sz(1);
        N=sz(2);
        MagEx=ones(M,N);
        PhaseEx=ones(M,N);
        NumPoint = jj;
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
        MagEx = mag2db(MagEx);
    % Uncomment to do figure analysis i.e. ellipsis
        MagEx = imbinarize(MagEx, -3);
        info = regionprops(MagEx, "Eccentricity","MajorAxisLength","MinorAxisLength");
        eccen(jj,ii) = info.Eccentricity;
        maj_ax(jj,ii) = info.MajorAxisLength;
        min_ax(jj,ii) = info.MinorAxisLength;
    end
    
end
%%
figure;
hold on;
for i=1:8
    plot([0,20,30,40,45,50],eccen(i,:),'--','Marker','*','LineWidth',1.5, ...
        'MarkerSize',8);
end
hold off;
grid on;
xlabel("angle (degrees)");
ylabel("eccentricity");
title("Eccentricity");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
legend(["250GHz","260GHz","270GHz","280GHz","290GHz","300GHz","310GHz","320GHz"],"Location","southeast");


figure;
hold on;
for i=1:8
    plot([0,20,30,40,45,50],maj_ax(i,:)*0.5,'--','Marker','*','LineWidth',1.5, ...
        'MarkerSize',8);
end
hold off;
grid on;
xlabel("angle (degrees)");
ylabel("major axis (mm)");
title("Major Axis");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
legend(["250GHz","260GHz","270GHz","280GHz","290GHz","300GHz","310GHz","320GHz"],"Location","northwest");


figure;
hold on;
for i=1:8
    plot([0,20,30,40,45,50],min_ax(i,:)*0.5,'--','Marker','*','LineWidth',1.5, ...
        'MarkerSize',8);
end
hold off;
grid on;
xlabel("angle (degrees)");
ylabel("minor axis (mm)");
title("Minor Axis");
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
legend(["250GHz","260GHz","270GHz","280GHz","290GHz","300GHz","310GHz","320GHz"],"Location","southeast");






