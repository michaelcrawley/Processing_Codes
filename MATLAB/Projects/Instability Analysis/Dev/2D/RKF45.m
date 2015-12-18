
function [z,y]=RKF45(z1,ystart,yend,abserr,alpha,c)
%   this code integrates a system of first order ordinary diferential
%   equations by runge-kutta-fehberg-45 method with automatic estimation
%   of local error and step size adjustment.
%   number of equations in the ODE system = size of z1.
%
%   [input]
%   z1 - vector of initial values of the dependent variables.
%   ystart - start value of the independent variable.
%   yend - end value of the independent variable.
%       When ystart>yend, use h= -0.05; while y(j)>yend; if y(j)+h<yend;
%       When ystart<yend, use h=  0.05; while y(j)<yend; if y(j)+h>yend;
%   abserr - bound of local error permitted, the computed solution
%       znew obtained in a step h must pass the test to be accepted
%       err<abserr
%   alpha - to be used in subfunction deriv, comes from the main program
%   c - to be used in subfunction deriv, comes from the main program
%
%   [output]
%   z - matrix of dependent variables at every integration point,
%       a column of z corresponds to a component in the y vector
%   y - vector of independent variable at every integration point.
%
%   [Comments]
%(1) The part of checking accuray and finding new h comes from document 10/16/2006
%(2) It's possible that: scheme for new h can not fulfill the error tolerance
%(3) To avoid dead loop, set maximun number of valid/total integrations.
%    When reaching maximun number of integration, go out of the function

% parameters of RKF45 method
a2=1/4; b2=1/4;
a3=3/8; b3=3/32; c3=9/32;
a4=12/13; b4=1932/2197; c4=-7200/2197; d4=7296/2197;
a5=1; b5=439/216; c5=-8; d5=3680/513; e5=-845/4104;
a6=1/2; b6=-8/27; c6=2; d6=-3544/2565; e6=1859/4104; f6=-11/40;
n1=25/216; n3=1408/2565; n4=2197/4104; n5=-1/5;
r1=1/360; r3=-128/4275; r4=-2197/75240; r5=1/50; r6=2/55;

% step size
h= 0.05;      % h can be either positive or negative
hmin=0.001;   % hmin must be positive
hmax=0.05;    % hmax must be positive

% control parameters to avoid dead loop
max1=1000;    % maximun number of valid integrations
max2=2000;    % maximun number of all integrations
i=0;          % number of all integration (including not valid integrations)
j=1;          % number of valid integrations

% results to be output
y=[ ];
y=[y,ystart];
z=[ ];
z=[z, z1];

while y(j)<yend
    if y(j)+h>yend
        h=yend-y(j);
    end
    k1=h*deriv(y(j),z(:,j), alpha,c);
    k2=h*deriv(y(j)+a2*h, z(:,j)+b2*k1, alpha,c);
    k3=h*deriv(y(j)+a3*h, z(:,j)+b3*k1+c3*k2, alpha, c);
    k4=h*deriv(y(j)+a4*h, z(:,j)+b4*k1+c4*k2+d4*k3, alpha, c);
    k5=h*deriv(y(j)+a5*h, z(:,j)+b5*k1+c5*k2+d5*k3+e5*k4, alpha, c);
    k6=h*deriv(y(j)+a6*h, z(:,j)+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5, alpha, c);
    % new z at new y point
    znew=z(:,j)+n1*k1+n3*k3+n4*k4+n5*k5;
    % local error
    err=abs(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6);
    err=max(err);
    % if err<abserr, accept new z; if not, go on to try new h
    if err<abserr
        y=[y,y(j)+h];
        z=[z,znew];
        j=j+1;
    end
    % new h (always s>0)
    s=(abserr*abs(h)/(2*err))^0.25;
    if s<0.1
        s=0.1;
    elseif s>4.0
        s=4.0;
    end
    h=s*h;
    if abs(h)>hmax
        h=hmax*h/abs(h);
    elseif abs(h)<hmin
        h=hmin*h/abs(h);
    end
    % reach maximum number of integration, break
    i=i+1;
    if j==max1 || i==max2
        disp(' ');
        disp('error in RKF45: too many integrations needed');
        return
    end
end






% Subfunctions
%=================================
function dz=deriv(y,z,alpha,c)
% y - independent variable
% z - dependent variable
% alpha - complex wave number
% c = omega/alpha
% dz = dz/dy at y
dz=zeros(4,1);
[u,u2]=u_sub(y);
cmc=u-c;
dz(1)=z(2);                        
dz(2)=(alpha*alpha+u2/cmc)*z(1) ;  
dz(3)=z(4);						   
dz(4)=(alpha*alpha+u2/cmc)*z(3)+z(1)*(2*alpha-c*u2/(cmc*cmc*alpha));




%=================================
function [u,u2]=u_sub(y)
% y - independent variable
% u - velocity at y
% u2 - second derivative of u at y
L=0.45;
W=0.8;
u=1+L*tanh(y)-W*exp(-log(2)*y.^2);
u2=-2*L*tanh(y).*(1-tanh(y).^2)+2*W*log(2)*exp(-log(2)*y.^2)-4*W*log(2)^2*y.^2.*exp(-log(2)*y.^2);


