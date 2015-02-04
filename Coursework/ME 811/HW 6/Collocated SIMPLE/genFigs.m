clear; close all; clc;

load re100.mat

x=linspace(0,L,N); y=linspace(0,H,M);

%% Plot the residuals
figure;
semilogy(R2.x,'-k'); hold on;
semilogy(R2.y,'--k');
semilogy(R2.p,'-k','LineWidth',1); hold off;
xlabel('Iterations'); ylabel('L_2 Norm of Residual');
legend('u','v','p');

%% Plot contours
figure;
% Plot the u-velocity
subplot(3,1,1); contourf(x,y,u.c);
colormap(flipud(gray)); c=colorbar;
axis equal; xlim([0 L]); ylim([0 H]);
ylabel('y [m]'); set(get(c,'ylabel'),'String','u [m/s]');

% Plot the v-velocity
subplot(3,1,2); contourf(x,y,v.c);
colormap(flipud(gray)); c=colorbar;
axis equal; xlim([0 L]); ylim([0 H]);
ylabel('y [m]'); set(get(c,'ylabel'),'String','v [m/s]');

% Plot the pressure
subplot(3,1,3); contourf(x,y,p.c);
colormap(flipud(gray)); c=colorbar;
axis equal; xlim([0 L]); ylim([0 H]);
xlabel('x [m]'); ylabel('y [m]');
set(get(c,'ylabel'),'String','p_{gauge} [Pa]');

%% Plot a vector graph
figure; quiver(x,y,u.c,v.c,2,'-k');
axis equal; xlim([0 L]); ylim([0 H]);
xlabel('x [m]'); ylabel('y [m]');