mainname = 'Fullstep_Extrange';

%2CDS1 
filename = [mainname,'_2CDS1'];
plot(x1,u2CDS1(end,:),x1,ua1,'-.');xlabel('x');ylabel('A');ylim([-0.2 0.5]);legend('2CDS','Analytic','Location','Northwest');
xlim([0 x1(end)]);
title('2^n^d Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x1,u2CDS1(end,:)-ua1);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x1(end)]);
title('Error in 2^n^d Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%4CDS1
filename = [mainname,'_4CDS1'];
plot(x1,u4CDS1(end,:),x1,ua1,'-.');xlabel('x');ylabel('A');ylim([-0.2 0.5]);legend('4CDS','Analytic','Location','Northwest');
xlim([0 x1(end)]);
title('4^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x1,u4CDS1(end,:)-ua1);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x1(end)]);
title('Error in 4^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%6CDS1
filename = [mainname,'_6CDS1'];
plot(x1,u6CDS1(end,:),x1,ua1,'-.');xlabel('x');ylabel('A');ylim([-0.2 0.5]);legend('6CDS','Analytic','Location','Northwest');
xlim([0 x1(end)]);
title('6^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x1,u6CDS1(end,:)-ua1);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x1(end)]);
title('Error in 6^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%4DRP1
filename = [mainname,'_4DRP1'];
plot(x1,u4DRP1(end,:),x1,ua1,'-.');xlabel('x');ylabel('A');ylim([-0.2 0.5]);legend('4DRP','Analytic','Location','Northwest');
xlim([0 x1(end)]);
title('4^t^h Order DRP Optimized Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x1,u4DRP1(end,:)-ua1);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x1(end)]);
title('Error in 4^t^h Order DRP Optimized Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%2CDS2 
filename = [mainname,'_2CDS2'];
plot(x2,u2CDS2(end,:),x2,ua2,'-.');xlabel('x');ylabel('A');ylim([-0.2 1.4]);legend('2CDS','Analytic','Location','Northwest');
xlim([0 x2(end)]);
title('2^n^d Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x2,u2CDS2(end,:)-ua2);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x2(end)]);
title('Error in 2^n^d Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%4CDS2
filename = [mainname,'_4CDS2'];
plot(x2,u4CDS2(end,:),x2,ua2,'-.');xlabel('x');ylabel('A');ylim([-0.2 1.4]);legend('4CDS','Analytic','Location','Northwest');
xlim([0 x2(end)]);
title('4^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x2,u4CDS2(end,:)-ua2);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x2(end)]);
title('Error in 4^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%6CDS2
filename = [mainname,'_6CDS2'];
plot(x2,u6CDS2(end,:),x2,ua2,'-.');xlabel('x');ylabel('A');ylim([-0.2 1.4]);legend('6CDS','Analytic','Location','Northwest');
xlim([0 x2(end)]);
title('6^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x2,u6CDS2(end,:)-ua2);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x2(end)]);
title('Error in 6^t^h Order Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

%4DRP2
filename = [mainname,'_4DRP2'];
plot(x2,u4DRP2(end,:),x2,ua2,'-.');xlabel('x');ylabel('A');ylim([-0.2 1.4]);legend('4DRP','Analytic','Location','Northwest');
xlim([0 x2(end)]);
title('4^t^h Order DRP Optimized Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');

filename = [filename,'e'];
plot(x2,u4DRP2(end,:)-ua2);xlabel('x');ylabel('error (Analytic-Numeric)');
xlim([0 x2(end)]);
title('Error in 4^t^h Order DRP Optimized Central Difference, \Deltax = 0.5, \Deltat = 0.05, Free boundaries')
saveas(gcf,filename,'fig');saveas(gcf,filename,'png');