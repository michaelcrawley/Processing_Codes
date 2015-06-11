function [Mimb R2] = calcMimb(u,v,rho)
% Calculate the mass imbalance.

M = size(u,1);
N = size(u,2)-1;
dy=0.01/M; dx=0.01/N;

i=1:M; j=1:N;
Mimb = rho*(u(i,j+1)-u(i,j))*dy + rho*(v(i+1,j)-v(i,j))*dx;

R2 = sqrt(sum(sum( Mimb.^2 )));