function [u v] = correctUV(u,v,dp,Aox,Aoy,w)

M = size(u,1);
N = size(u,2)-1;
dy=0.01/M; dx=0.01/N;

%% Initialize matrices
du = zeros(M,N+1);
dv = zeros(M+1,N);

%% Calculate corrections
i=1:M; j=2:N;
du(i,j) = ( dp(i,j-1)-dp(i,j) )*dy./Aox(i,j);

i=2:M; j=1:N;
dv(i,j) = ( dp(i-1,j)-dp(i,j) )*dx./Aoy(i,j);

%% Apply corrections with relaxation factor
u = u + w*du;
v = v + w*dv;