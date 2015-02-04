load temp_oldBLdata.mat
Po=48; %psig
rem=6; %number of points to remove from the bottom of the profile
u=u(rem:end);
y=y(rem:end);
utau=utau*3.125;

%function [yn,un]=mmBLplot(y,delta,Minf,u,uinf,utau,Po)

%Calculates and plots the van Driest transformed, normalized velocity
%profile as proposed by Maise and McDonald in their paper "Mixing Length
%and Kinematic Eddy Viscosity in a Compressible Boundary Layer." 1968
%
%The function input/output is as follows:
%
%[yn,un]=mmBLplot(y,delta,Minf,u,uinf,utau)
%
%where yn and un are the coordinate and velocity of the normalized profile,
%y is the dimensional coordinate, delta is the boundary layer thickness
%(must be dimensionally consisten with y), Minf is the freestream Mach
%number, u is the dimensional velocity profile, uinf is the freestream
%velocity, and utau is the friction velocity. Note that all velocities MUST
%be in m/s.

Tinf=(uinf/Minf)^2/(1.4*287);
To=Tinf*(1+0.2*Minf^2);
Tw=(To+Tinf)/2;


% %EXPERIMENTAL SECTION
% %calculating the friction velocity
% N=7; %index of last point over which to average dudy at the wall
% 
% %viscosity is calculated using southerland's law
% Tos=524.07; %degR
% muos=0.01827; %centipoise
% C=120; %southerland's constant
% TR=((9/5*(Tw+273.15))+32)-459.67;
% a=0.555*Tos+C;
% b=0.555*TR+C;
% mu=0.001*muos*(a/b)*(TR/Tos)^(3/2);
% 
% %density calculation
% Po=101325*(Po+14.696)/14.696;
% Pinf=Po*(1+0.2*Minf^2)^(-3.5);
% rho=Pinf/(287*Tw);
% 
% dudy=0;
% for i=2:N
%     dudy=dudy+(u(i)-u(i-1))/((y(i)-y(i-1))/1000);
% end
% dudy=dudy/(N-1);
% 
% utau=sqrt(mu*dudy/rho);
% 
% %EXPERIMENTAL SECTION


%model profile
ymod=[0.01:0.01:1];
umod=-2.5.*log(ymod)+1.25.*(1+cos(3.14159265.*ymod));


%normalizing the experimental profile
yn=y./delta;
yn=yn(find(yn<=1));
u=u(1:length(yn));
A=sqrt(0.2*Minf^2/(Tw/Tinf));
B=(1+0.2*Minf^2)/(Tw/Tinf)-1;
us=uinf/A.*asin((2*A^2.*(u./uinf)-B)/sqrt(B^2+4*A^2));
uinfs=uinf/A.*asin((2*A^2-B)/sqrt(B^2+4*A^2));
un=(uinfs-us)./utau;


%calculating the error between the actual and model profiles
umodinterp=interp1(ymod,umod,yn,'spline');
err=(un-umodinterp)./umodinterp;


%Un-transforming the model profile
umods=uinfs-umod.*utau;
umodu=B/(2*A^2)+sqrt((B/(2*A^2))^2+A^-2).*sin(A.*umods./uinf);


% %trying to fit a utau to the profile by least squares regression
% utau_old=utau;
% low=1;
% high=length(un);
% llimutau=0.2;
% ulimutau=100;
% Ns=100;
% range=[llimutau:(ulimutau-llimutau)/(Ns-1):ulimutau];
% for i=1:Ns
%     d
% end


%plotting the van Driest transformed normalized profiles
figure
set(gcf,'Color','White');
semilogx(ymod,umod,'-k',yn,un,'+k');
title('van Driest Transformed, Normalized Velocity Profile');
xlabel('y/\delta');
ylabel('(u^*_\infty-u^*)/u_\tau');
legend('Model','Actual');
set(legend,'Location','NorthEast');

%plotting the normalized profiles
figure
set(gcf,'Color','White');
plot(umodu,ymod,'-k',u./uinf,yn,'+k');
set(gca,'XLim',[0,1]);
set(gca,'YLim',[0,1]);
title('Normalized Velocity Profiles')
xlabel('u/u_\infty');
ylabel('y/\delta');
legend('Model','Actual');
set(legend,'Location','NorthEast');

%plotting the error in the van Driest transformed velocity profiles
figure
set(gcf,'Color','White');
semilogx(yn,err,'+k');
set(gca,'YLim',[-0.4,0.1]);
set(gca,'XLim',[0.1,1]);
xlabel('Normalized Y Coordinate');
ylabel('Relative Error in the Experimental Profile');
title('Error in the Experimental Profile Relative to the Model');