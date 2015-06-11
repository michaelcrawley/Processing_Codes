function [Mimb R2] = calcMimb(uf,vf,rhof,dy,dx)
% Calculate the mass imbalance.

M = size(uf,1);
N = size(vf,2);

% Initialize matrix
Mimb = zeros(M,N);

% Calculate mass-imbalance for all cells
i=1:M; j=1:N;
Mimb = ( rhof.x(i,j+1).*uf(i,j+1) - rhof.x(i,j).*uf(i,j) )*dy ...
	+ ( rhof.y(i+1,j).*vf(i+1,j) - rhof.y(i,j).*vf(i,j) )*dx;

% Calculate the residual
R2 = sqrt(sum(sum( Mimb.^2 )));