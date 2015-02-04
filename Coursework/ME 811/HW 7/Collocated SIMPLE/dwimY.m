function [fy] = dwimY(c,M,N,m,n)
% Calculates the distance weighted interpolation of cell face values
% assuming a uniform grid. Domain boundaries are assumed equal to nearest
% node value. Step boundaries are not handled.

fy = zeros(M+1,N);

% Bottom boundaries
i=1; j=n+1:N;
fy(i,j) = c(i,j);
i=m+1; j=1:n;
fy(i,j) = c(i,j);

% Interior boundaries
i=2:m+1; j=n+1:N;
fy(i,j) = (c(i-1,j) + c(i,j))/2;
i=m+2:M; j=1:N;
fy(i,j) = (c(i-1,j) + c(i,j))/2;

% Top boundary
i=M+1; j=1:N;
fy(i,j) = c(i-1,j);