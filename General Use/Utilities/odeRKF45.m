function [z,y]=odeRKF45(fun,z1,ystart,yend,abserr,inputs)
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

if yend < ystart
    dx = -1;
else
    dx = 1;
end
% step size
h= 0.01*dx;      % h can be either positive or negative
hmin=0.0001*dx;   % hmin must be positive
hmax=0.02*dx;    % hmax must be positive

% control parameters to avoid dead loop
max1=200000;    % maximun number of valid integrations
max2=400000;    % maximun number of all integrations
i=0;          % number of all integration (including not valid integrations)
j=1;          % number of valid integrations

% results to be output
y=[ ];
y=[y,ystart];
z=[ ];
z=[z, z1];

while abs(y(j)-yend) > eps
    if abs(y(j)-yend) < abs(h)
        h=abs(yend-y(j))*dx;
    end
    k1=h*feval(fun,y(j),z(:,j), inputs);
    k2=h*feval(fun,y(j)+a2*h, z(:,j)+b2*k1, inputs);
    k3=h*feval(fun,y(j)+a3*h, z(:,j)+b3*k1+c3*k2, inputs);
    k4=h*feval(fun,y(j)+a4*h, z(:,j)+b4*k1+c4*k2+d4*k3, inputs);
    k5=h*feval(fun,y(j)+a5*h, z(:,j)+b5*k1+c5*k2+d5*k3+e5*k4, inputs);
    k6=h*feval(fun,y(j)+a6*h, z(:,j)+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5, inputs);
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
        h=dx*hmax*h/abs(h);
    elseif abs(h)<hmin
        h=dx*hmin*h/abs(h);
    end
    % reach maximum number of integration, break
    i=i+1;
    if j==max1 || i==max2
        %disp(' ');
        %disp('error in RKF45: too many integrations needed');
        return
    end

end