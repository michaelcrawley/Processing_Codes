function [fx] = dwimX(c,M,N,m,n)
% Calculates the distance weighted interpolation of cell face values
% assuming a uniform grid. Domain boundaries are assumed equal to nearest
% node value. Step boundaries are not handled.

fx = zeros(M,N+1);

% Left boundaries
i=m+1:M; j=1;
fx(i,j) = c(i,j);
i=1:m; j=n+1;
fx(i,j) = c(i,j);

% Interior boundaries
i=m+1:M; j=2:n+1;
fx(i,j) = (c(i,j-1) + c(i,j))/2;
i=1:M; j=n+2:N;
fx(i,j) = (c(i,j-1) + c(i,j))/2;

% Right boundary
i=1:M; j=N+1;
fx(i,j) = c(i,j-1);