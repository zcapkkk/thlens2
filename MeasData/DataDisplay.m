%Calculate the radiation pattern 
clear all;close all;
clc;


fname = '20211231_B50deg.mat';
load(fname);
for NumPoint = 1:7
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
    title_string = {'Amp (a.u),',int2str(240+10*NumPoint),'GHz,',fname(11:12),'Deg.'};
    title_string = strjoin(title_string);
    file_string = {'./Results/','freq',int2str(240+10*NumPoint),'GHz',fname(11:12),'.png'};
    file_string = strjoin(file_string,'');
    title(title_string);
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
    saveas(gcf, file_string);
end

