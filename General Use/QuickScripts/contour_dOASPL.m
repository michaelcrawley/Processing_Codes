R = (-5:0.5:10);
ChPol = [90 80 70 60 50 45 40 35 30 25];

CX = [-4 8];
TTR = '2.5';
%%%% TOTAL SPECTRUM %%%%
figure
[C,H] = contourf(ChPol,F,dO(:,:,1),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL (dB), m = 0, T_o/T_a = ' TTR])
saveas(gcf,'dO_m0.fig')
saveFigure_v2(gcf,'dO_m0',600)

figure
[C,H] = contourf(ChPol,F,dO(:,:,2),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL (dB), m = 1, T_o/T_a = ' TTR])
saveas(gcf,'dO_m1.fig')
saveFigure_v2(gcf,'dO_m1',600)

figure
[C,H] = contourf(ChPol,F,dO(:,:,3),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL (dB), m = 3, T_o/T_a = ' TTR])
saveas(gcf,'dO_m3.fig')
saveFigure_v2(gcf,'dO_m3',600)

figure
[C,H] = contourf(ChPol,F,dO(:,:,4),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL (dB), m = \pm4, T_o/T_a = ' TTR])
saveas(gcf,'dO_m4.fig')
saveFigure_v2(gcf,'dO_m4',600)




%%%% DETONED SPECTRUM %%%%
R = (0:0.3:5); R = [fliplr(-R(2:end)) R];

CX = [-5 5];
figure
[C,H] = contourf(ChPol,F,dOT(:,:,1),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL - Detoned (dB), m = 0, T_o/T_a = ' TTR])
saveas(gcf,'dOT_m0.fig')
saveFigure_v2(gcf,'dOT_m0',600)

figure
[C,H] = contourf(ChPol,F,dOT(:,:,2),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL - Detoned (dB), m = 1, T_o/T_a = ' TTR])
saveas(gcf,'dOT_m1.fig')
saveFigure_v2(gcf,'dOT_m1',600)

figure
[C,H] = contourf(ChPol,F,dOT(:,:,3),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL - Detoned (dB), m = 3, T_o/T_a = ' TTR])
saveas(gcf,'dOT_m3.fig')
saveFigure_v2(gcf,'dOT_m3',600)

figure
[C,H] = contourf(ChPol,F,dOT(:,:,4),R);
xlabel('Polar Angle (degrees)'); set(gca,'XDir','reverse');
ylabel('St_{DF}');
grid on
axis([25 90 0.09 3])
colorbar
caxis(CX);
colormap(CM)
clabel(C,H,'LabelSpacing',72*3,'Rotation',0,'FontWeight','bold','Color',[255 100 255]/255);
title(['\DeltaOASPL - Detoned (dB), m = \pm4, T_o/T_a = ' TTR])
saveas(gcf,'dOT_m4.fig')
saveFigure_v2(gcf,'dOT_m4',600)
