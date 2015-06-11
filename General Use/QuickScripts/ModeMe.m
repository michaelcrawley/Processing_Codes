clear variables

dir_name = 'F:\Research\RingGrooveEffect\20091209\';
O1 = 5; %offset 1
O2 = 15;    %offset 2
YL = [40 85];   %y-axis limits
TTR = 1.0;  %total temperature ratio
NP = 29;    %number of forcing frequencies

flist = struct2cell(dir(dir_name));
DM = flist(2,:);
q = strfind(flist(1,:),'S.fftNOS'); %extracts list of data files - ignores all other files
keep = ~cellfun('isempty',q);
flist = flist(1,logical(keep))';
BL = flist(~cellfun('isempty',strfind(flist,'Baseline')));
mL{1} = [];
mL(1:NP,1) = flist(~cellfun('isempty',strfind(flist,'m(00)')));
mL(1:NP,2) = flist(~cellfun('isempty',strfind(flist,'m(01)(01)')));
mL(1:NP,3) = flist(~cellfun('isempty',strfind(flist,'m(02)(02)')));
mL(1:NP,4) = flist(~cellfun('isempty',strfind(flist,'m(03)')));
mL(1:NP,5) = flist(~cellfun('isempty',strfind(flist,'m(01)(-1)')));
mL(1:NP,6) = flist(~cellfun('isempty',strfind(flist,'m(02)(-2)')));
mL(1:NP,7) = flist(~cellfun('isempty',strfind(flist,'m(04)(-4)')));

for n = 1:NP
    [Names,Vals,kp] = readSETfile([dir_name mL{n,1}(1:end-9) '.SET'],[38 39]); 
    F(n,1) = Vals{1}/1000; STdf(n,1) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,2}(1:end-9) '.SET'],[38 39]); 
    F(n,2) = Vals{1}/1000; STdf(n,2) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,3}(1:end-9) '.SET'],[38 39]); 
    F(n,3) = Vals{1}/1000; STdf(n,3) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,4}(1:end-9) '.SET'],[38 39]); 
    F(n,4) = Vals{1}/1000; STdf(n,4) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,5}(1:end-9) '.SET'],[38 39]); 
    F(n,5) = Vals{1}/1000; STdf(n,5) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,6}(1:end-9) '.SET'],[38 39]); 
    F(n,6) = Vals{1}/1000; STdf(n,6) = Vals{2};
    [Names,Vals,kp] = readSETfile([dir_name mL{n,7}(1:end-9) '.SET'],[38 39]); 
    F(n,7) = Vals{1}/1000; STdf(n,7) = Vals{2};
end
for n = 1:size(mL,2)
    [STdf(:,n),IX] = sort(STdf(:,n)); F(:,n) = F(IX,n); mL(:,n) = mL(IX,n);
end
% STdf = [STdf(:,1) STdf(:,2) STdf(:,4) STdf(:,7)];
STdf = mean(STdf,2);
% F = [F(:,1) F(:,2) F(:,4) F(:,7)];
F = mean(F,2);

B = zeros(4086,12);
for n = 1:length(BL)
    B = B + dlmread([dir_name BL{n}],'\t', 1, 0);
end
B = B/n;

for n = 1:NP
    disp(n)
    m0 = dlmread([dir_name mL{n,1}],'\t', 1, 0);
    m1 = dlmread([dir_name mL{n,2}],'\t', 1, 0);
    m2 = dlmread([dir_name mL{n,3}],'\t', 1, 0);
    m3 = dlmread([dir_name mL{n,4}],'\t', 1, 0);
    m4 = dlmread([dir_name mL{n,5}],'\t', 1, 0);
    m5 = dlmread([dir_name mL{n,6}],'\t', 1, 0);
    m6 = dlmread([dir_name mL{n,7}],'\t', 1, 0);
    
    figure
    semilogx(B(:,2),B(:,3),'-.k','LineWidth',1);
    hold on
    semilogx(m0(:,2),m0(:,3),'-r','LineWidth',1);
    semilogx(m1(:,2),m1(:,3),'-g','LineWidth',1);
    semilogx(m2(:,2),m2(:,3),'-b','LineWidth',1);
    semilogx(m3(:,2),m3(:,3),'-c','LineWidth',1);
    semilogx(m4(:,2),m4(:,3),'-m','LineWidth',1);
    semilogx(m5(:,2),m5(:,3),'-y','LineWidth',1);
    semilogx(m6(:,2),m6(:,3),'-k','LineWidth',1);
    
    semilogx(B(:,2),B(:,7)+O1,'-.k','LineWidth',1);
    semilogx(m0(:,2),m0(:,7)+O1,'-r','LineWidth',1);
    semilogx(m1(:,2),m1(:,7)+O1,'-g','LineWidth',1);
    semilogx(m2(:,2),m2(:,7)+O1,'-b','LineWidth',1);
    semilogx(m3(:,2),m3(:,7)+O1,'-c','LineWidth',1);
    semilogx(m4(:,2),m4(:,7)+O1,'-m','LineWidth',1);
    semilogx(m5(:,2),m5(:,7)+O1,'-y','LineWidth',1);
    semilogx(m6(:,2),m6(:,7)+O1,'-k','LineWidth',1);

    
    semilogx(B(:,2),B(:,11)+O2,'-.k','LineWidth',1);
    semilogx(m0(:,2),m0(:,11)+O2,'-r','LineWidth',1);
    semilogx(m1(:,2),m1(:,11)+O2,'-g','LineWidth',1);
    semilogx(m2(:,2),m2(:,11)+O2,'-b','LineWidth',1);
    semilogx(m3(:,2),m3(:,11)+O2,'-c','LineWidth',1);
    semilogx(m4(:,2),m4(:,11)+O2,'-m','LineWidth',1);
    semilogx(m5(:,2),m5(:,11)+O2,'-y','LineWidth',1);
    semilogx(m6(:,2),m6(:,11)+O2,'-k','LineWidth',1);

    hold off
    legend('Base','m = 0','m = 1','m = 2','m = 3','m = \pm1','m = \pm2','m = \pm4','Location','NorthWest')

    axis([0.01 4 YL])
    grid on
    xlabel('St_D')
    ylabel('SPL (dB)')
    set(gca,'yTickLabel',[])
%     YL = get(gca,'YLim');
    set(gca,'YTick',(YL(1):5:YL(2)))
    yTickSeparationMarker(0.001,5,0.01,50)

    title({['T_o/T_a = ' num2str(TTR,'%.1f') ', F = ' num2str(F(n),'%.2f') ' kHz, St_{DF} = ' num2str(STdf(n),'%.2f')],['90^o\rightarrow+0, ' ...
        '50^o\rightarrow+' num2str(O1) ', 30^o\rightarrow+' num2str(O2)]})


    saveas(gcf,[cd '\TTR' num2str(TTR,'%.1f') '_STdf' num2str(STdf(n),'%.2f') '.fig'],'fig')
    saveFigure_v2(gcf,[cd '\TTR' num2str(TTR,'%.1f') '_STdf' num2str(STdf(n),'%.2f')],300)
    close
end
