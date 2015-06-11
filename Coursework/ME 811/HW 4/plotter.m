filename = [num2str(M),'x',num2str(N),' mesh, Re = ',num2str(round(rho*Ulid*L/gamma))];

%Plot centerline u,v velocity along centerlines
h(1) = figure; %y = const figure
plot(x,(u(M/2,:)+u(M/2+1,:))/2,'k',xc,v(M/2+1,:),'k--');legend('u velocity','v velocity','Location','Best');xlabel('x');ylabel('(m/s)');title(['Horizontal Centerline Velocities for ',filename]);
H(2) = figure; %x = const figure
plot(yc,u(:,N/2+1),'k',y,(v(:,N/2)+v(:,N/2+1))/2,'k--');legend('u velocity','v velocity','Location','Best');xlabel('y');ylabel('(m/s)');title(['Vertical Centerline Velocities for ',filename]);

%Plot residuals
h(3) = figure;
semilogy(1:counter-1,Resx,'k',1:counter-1,Resy,'k--',1:counter-1,Resp,'k-.');legend('u Residual','v Residual','p Residual');xlabel('iteration');ylabel('Residual');title(['Outer Iteration Residuals for ',filename]);


%Plot quiver and pressure
uc = (u(:,1:end-1)+u(:,2:end))/2;
vc = (v(1:end-1,:)+v(2:end,:))/2;
h(4) = figure;
contourf(xc,yc,p,50);colormap(flipud(gray));colorbar;hold on;
quiver(xc,yc,uc,vc,3);xlabel('x');ylabel('y');hold off;
title(['Velocity Vectors and Pressure Contours for ',filename]);

%Save figures
saveas(h(1),[filename,' y centerline'],'fig');saveas(h(1),[filename,' y centerline'],'png');
saveas(h(2),[filename,' x centerline'],'fig');saveas(h(2),[filename,' x centerline'],'png');
saveas(h(3),[filename,' Residuals'],'fig');saveas(h(3),[filename,' Residuals'],'png');
saveas(h(4),[filename,' pressure'],'fig');saveas(h(4),[filename,' pressure'],'png');

close(h);