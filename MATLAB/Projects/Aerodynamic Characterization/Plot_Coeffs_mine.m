function Plot_Coeffs(aero,color)

%figure();
subplot(3,3,1);
plot([aero.AoA], [aero.CL],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Lift');
subplot(3,3,2);
plot([aero.AoA], [aero.CD],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Drag');
subplot(3,3,3);
plot([aero.AoA], [aero.CS],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Side');
subplot(3,3,4);
plot([aero.AoA], [aero.CRM],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Roll');
subplot(3,3,5);
plot([aero.AoA], [aero.CPM],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Pitch');
subplot(3,3,6);
plot([aero.AoA], [aero.CYM],color);
hold on;
grid on;
xlabel('AoA');
ylabel('C Yaw');

subplot(3,3,7);
plot([aero.AoA], [aero.CL]./[aero.CD],color);
hold on;
grid on;
xlabel('AoA');
ylabel('L/D');

subplot(3,3,8);
plot([aero.AoA], [aero.CL].^1.5./[aero.CD],color);
hold on;
grid on;
xlabel('AoA');
ylabel('CL^{1.5}/CD');

subplot(3,3,9);
plot([aero.CD], [aero.CL],color);
hold on;
grid on;
xlabel('CD');
ylabel('CL');

saveas(gcf,'All_Data.fig');
export_fig('All_Data','-png','-r300');

