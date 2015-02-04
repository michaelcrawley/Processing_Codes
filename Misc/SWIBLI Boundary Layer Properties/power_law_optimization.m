%temporary back-of-the-envelope n calculation
close all
u_n=[.4078,.6531,.8449,.9264,.9607,.981,.9901];%
y=[.4812,1.283,2.085,2.887,3.689,4.491,5.293];%

u=u_n.*482; %m/s
y_n=y./4.8; %mm

n=bl_power_law_fit(y,u);

ymodn=[0:0.01:5]./4.8;
umodn=ymodn.^(1/n);

%calculating 1/7 law profile
y7n=ymodn;
u7=y7n.^(1/7);
fprintf('\n\nBest fit curve is with n = %f\n',n);

%plotting the boundary layer profiles
figure
set(gcf,'Color','White');
plot(u_n,y_n,'+k',umodn,ymodn,'-k',u7,y7n,'--k');
title('Experimental, Best-Fit and 1/7 Power Law Profiles');
xlabel('u/u_\infty');
set(gca,'XLim',[0,1]);
ylabel('y/\delta');
set(gca,'YLim',[0,1.2]);
legend('Experimental','Best-Fit','1/7 Power Law');
set(legend,'Location','NorthWest');

%calculation incompressible shape factor

%calculation compressible shape factor
[cdist,dr]=cmp_disc_thickness(u,y./1000,482,1.89);
[cmomt,dr]=cmp_mom_thickness(u,y./1000,482,1.89);
Hcmp=dist/momt;
fprintf('The compressible displacement thickness is d* = %.3f mm\n',cdist*1000)
fprintf('The compressible momentum thickness is O = %.3f mm\n',cmomt*1000)
fprintf('The compressible shape factor is H = %.3f\n\n',Hcmp);

%plotting the density ratio profile
figure
set(gcf,'Color','White');
plot(dr,y_n,'-k');
xlabel('\rho/\rho_\infty');
ylabel('y/\delta');