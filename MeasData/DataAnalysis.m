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

for ii = 1:length(degrees)
    for jj = 1:9
        fname = convertStringsToChars(degrees(ii));
        load(fname);

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
%         MagEx=MagEx/max(MagEx(:));
%         MagEx_dB=20*log10(MagEx);
        MagEx_dB = 20*log10(MagEx);
% Uncomment to do figure analysis i.e. ellipsis
%         MagEx = imbinarize(MagEx, -5);
%         info = regionprops(MagEx, "Eccentricity","MajorAxisLength","MinorAxisLength");
%         diary res.out;
%         disp("'Deg': "+fname(11:12)+",")
%         disp("'Freq': "'+" "+int2str(240+10*NumPoint)+",");
%         disp("'Eccentricity': "+info.Eccentricity+",");
%         disp("'Maj. Axis': "+info.MajorAxisLength+",");
%         disp("'Min. Axis': "+info.MinorAxisLength+",");
%         diary off;
        %%%%%%%%%%%%%%%%%%%%%Figure Display%%%%%%%%%%%%%%%%%%%%
        figure('WindowStyle', 'Docked');
        x_a = -8:0.5:8;
        y_a = 50:-0.5:-10;
%         [max_r, max_c] = find(MagEx_dB==0);
        [max_r, max_c] = find(max(max(MagEx_dB)));
        y_a = y_a - y_a(max_r);
        imagesc(x_a,y_a,MagEx_dB);
        colorbar;
        colormap(jet);
        axis equal;
        xlabel('x (mm)');
        ylabel('y (mm)');
        axis xy;
        title_string = {'Amp (dB),',int2str(240+10*NumPoint),'GHz,',fname(11:12),'Deg.'};
        title_string = strjoin(title_string);
        file_string = {'./Results/','freq',int2str(240+10*NumPoint),'GHz',fname(11:12),'.png'};
        file_string = strjoin(file_string,'');
        title(title_string);
%         caxis([-20 0]);
        caxis([-90 0]);
        set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
%         saveas(gcf, file_string);
    end
end




