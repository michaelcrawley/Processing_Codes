function [u,u1]=u_sub(y)
% y - independent variable
% u - velocity at y
% u2 - second derivative of u at y
a=0.25*10;%6.25;
u=0.5*(1.+tanh(a*(1./y-y)));
u1=0.5*(1-(tanh(a*(1./y-y)))^2)*a*(-1./y^2-1.);
