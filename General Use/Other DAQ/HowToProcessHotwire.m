%How to process shear layer data

D = 0.0254*1.5; %m - jet diameter - remember to change as necessary in hotwireData.

    %Gather all the data
OFF = hotwireData('PSIG'); %special acquisition for establishing the origin of the coordinate system - assumes all data acquired on one coordinate system
M25 = hotwireData('PSIG');
M40 = hotwireData('PSIG');
M55 = hotwireData('PSIG');
M70 = hotwireData('PSIG');

    %Locate jet centerline
x = fliplr(linspace(min(OFF.POS),max(OFF.POS),length(OFF.POS)*20));
P = interp1(OFF.POS,OFF.MV,x);
Q = (min(P) + max(P))/2;
L = min(find(P >= Q));
R = max(find(P >= Q));
JC = (x(L)+x(R))/2;
clear x P Q L R

x = (JC-M40.POS)/1e6;   %create coordinates with proper zero - assumes DAQ positions are in microns

MV = [M25.MV' M40.MV' M55.MV' M70.MV']; %Gather profile voltages
load Calibrate.mat %load in hotwire calibration coefficients
for n = 1:4 %m/s - Convert voltage to velocity
    MV(:,n) = polyval(Cp,MV(:,n));
end
MV = MV./repmat(mean(MV(end-9:end,:)),154,1);    %Scale velocity by free stream - remember to change length of matrix.
clear n
save LAFPA.mat  %save data to mat file

    %Plot profiles
CP{1} = '-ob'; CP{2} = '-sr'; CP{3} = '-dg'; CP{4} = '-<k';
for N = 1:4
    plot(x/D,MV(:,N),CP{N},'LineWidth',1.5)
    hold on
end
grid on
xlabel('r/D')
ylabel('u/U_j')
axis([0.46 0.54 0 1])
legend('M = 0.25','M = 0.40','M = 0.55','M = 0.70')
title('LAFPA Extension Nozzle')
saveFigure(gcf,'LAFPA',600);    %remember to save as figure file as well

for N = 1:4 %Calculate shape parameters
    [yn(:,N),I(:,N),BL(:,N),DT(:,N),MT(:,N),H(:,N)] = BLcalc(x,MV(:,N));
end
clear N yn
%Remember to save this to mat file as well

    %Calculates momentum thickness using curve fitting
for n = 1:9
    [params(n,:), fCurve(:,n), sse(n)] = tanhfit(x, NV(:,n));
end
th = -1./params(:,4)/4;