%Calculate the radiation pattern 
close all;
clc;

deg00 = "20211230_B00deg.mat";
deg20 = "20220102_B20deg.mat";
deg30 = "20211231_B30deg.mat";
deg40 = "20211230_B40deg.mat";
deg45 = "20220102_B45deg.mat";
deg50 = "20211231_B50deg.mat";

degrees = [deg00, deg20, deg30, deg40, deg50];
store_oned = zeros(7,6,121);

for ii = 1:length(degrees)
    for jj = 1:8
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
        MagEx=MagEx/max(max(MagEx(:)));
        MagEx_dB=20*log10(MagEx);
        [maxrow, maxcol] = find(MagEx_dB==0);
        csec_MagEx_dB = MagEx_dB(:, maxcol);
        store_oned(jj,ii,:) = csec_MagEx_dB;    
    end
end


for i = 1:8
    figure;
    title_string = {'Amplitude Cross Section',int2str(240+10*i),'GHz'};
    title_string = strjoin(title_string);
    sv_str = strjoin({'./Amplitude_1D/1D_',int2str(240+10*i),'.png'},'');
    hold on;
    x_a = -10:0.5:50;
    for j = 1:5
        if j == 1
            zero_peak = x_a(reshape(flip(store_oned(i,j,:)),1,121)==0);
            x_a = x_a-zero_peak;
        end
        plot(x_a,reshape(flip(store_oned(i,j,:)),1,121),'LineWidth',2);
    end
    hold off;
    xlim([x_a(1) x_a(end)]);
    ylim([-20 0]);
    legend(["0$^\circ$","20$^\circ$","30$^\circ$","40$^\circ$","50$^\circ$"],'Interpreter','latex');
    xlabel("y (mm)");
    ylabel("dB");
%     title(title_string);
    grid MINOR;
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
    saveas(gcf,sv_str);
end


