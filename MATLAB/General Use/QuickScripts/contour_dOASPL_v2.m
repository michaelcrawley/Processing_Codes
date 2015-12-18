plot(F,m0.dAAE,'-or','LineWidth',1.5);
hold on
plot(F,m1.dAAE,'-^g','LineWidth',1.5);
plot(F,m2.dAAE,'-db','LineWidth',1.5);
plot(F,m3.dAAE,'-vc','LineWidth',1.5);
plot(F,Vm11.dAAE,'-sm','LineWidth',1.5);
plot(F,m22.dAAE,'-xy','LineWidth',1.5);
plot(F,m44.dAAE,'-*k','LineWidth',1.5);
hold off
grid on
xlabel('St_{DF}')
ylabel('\DeltaAAE (dB)')
legend('m = 0','m = 1','m = 2','m = 3','m = \pm1','m = \pm2','m = \pm4');
title('\DeltaAAE')

figure
plot(F,m0.dAAET,'-or','LineWidth',1.5);
hold on
plot(F,m1.dAAET,'-^g','LineWidth',1.5);
plot(F,m2.dAAET,'-db','LineWidth',1.5);
plot(F,m3.dAAET,'-vc','LineWidth',1.5);
plot(F,Vm11.dAAET,'-sm','LineWidth',1.5);
plot(F,m22.dAAET,'-xy','LineWidth',1.5);
plot(F,m44.dAAET,'-*k','LineWidth',1.5);
hold off
grid on
xlabel('St_{DF}')
ylabel('\DeltaAAE (dB)')
legend('m = 0','m = 1','m = 2','m = 3','m = \pm1','m = \pm2','m = \pm4');
title('\DeltaAAE - Detoned')




R = (-5:0.5:10);
TTR = '1.0';
%%%% TOTAL SPECTRUM %%%%
figure
[C,H] = contourf(ChPol,F,m0.dO,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL (dB), m = 0, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,dO(:,:,2),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL (dB), m = 1, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,dO(:,:,3),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL (dB), m = 3, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,dO(:,:,4),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL (dB), m = \pm4, T_o/T_a = ' TTR])





%%%% DETONED SPECTRUM %%%%
R = (-5:0.3:5);

figure
[C,H] = contourf(ChPol,F,m0.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = 0, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,m1.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = 1, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,m2.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = 2, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,m3.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = 3, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,Vm11.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = \pm1, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,Vm22.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = \pm2, T_o/T_a = ' TTR])

figure
[C,H] = contourf(ChPol,F,Vm44.dOT,R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color','m');
title(['\DeltaOASPL - Detoned (dB), m = \pm4, T_o/T_a = ' TTR])