%Calculate the radiation pattern 
clear all;close all;
clc;

NumPoint=7;
load('20220102_B20deg.mat');
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
%%%%%%%%%%%%%%%%%%%%%Figure Display%%%%%%%%%%%%%%%%%%%%
figure;
imagesc(MagEx);
colorbar;
colormap(jet);
axis equal;
xlabel('x');
ylabel('y');
title('Amp (a.u.)');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);

