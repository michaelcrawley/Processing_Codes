clear all; close all; clc;

load n40re100.mat

% Vertical centerline
figure;
plot([0 0],[0 0.01],'Color',[0.8 0.8 0.8],'HandleVisibility','off'); hold on;
plot(uc(:,21),dy/2:dy:0.01-dy/2,'-k');
plot(v(:,21),0:dy:0.01,'--k'); hold off;
xlabel('Velocity [m/s]'); ylabel('y [m]')
legend('u','v','Location','SouthEast');

% Horizontal centerline
figure;
plot([0 0.01],[0 0],'Color',[0.8 0.8 0.8],'HandleVisibility','off'); hold on;
plot(0:dx:0.01,u(21,:),'-k');
plot(dx/2:dx:0.01-dx/2,vc(21,:),'--k'); hold off;
xlabel('x [m]'); ylabel('Velocity [m/s]');
legend('u','v');

% Residual convergence
figure;
semilogy(R2.x,'-k'); hold on;
semilogy(R2.y,'--k');
semilogy(R2.p,'-k','LineWidth',1); hold off;
ylim(10.^[-10 -2]);
xlabel('Iterations'); ylabel('L_2-Normed Residual');
legend('u','v','p''');

% Vector plot (extra credit)
figure;
contourf(x,y,p); colormap(flipud(gray)); c = colorbar;
hold on; quiver(x,y,uc,vc,2,'-k'); hold off;
axis square; xlim([0 0.01]); ylim([0 0.01]);
xlabel('x [m]'); ylabel('y [m]');
set(get(c,'ylabel'),'String','Pressure (gauge), p [Pa]');