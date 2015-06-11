function [u v] = correctUV(u,v,dp,Ax,Ay,al,w)

M = size(u,1);
N = size(u,2)-1;
dy=0.01/M; dx=0.01/N;

%% Initialize matrices
du = zeros(M,N+1);
dv = zeros(M+1,N);

%% Calculate corrections
i=1:M; j=2:N;
du(i,j) = ( dp(i,j-1)-dp(i,j) )*dy ./ ((1+al)*Ax.o(i,j) + Ax.n(i,j) + Ax.s(i,j) + Ax.e(i,j) + Ax.w(i,j));

i=2:M; j=1:N;
dv(i,j) = ( dp(i-1,j)-dp(i,j) )*dx ./ ((1+al)*Ay.o(i,j) + Ay.n(i,j) + Ay.s(i,j) + Ay.e(i,j) + Ay.w(i,j));

%% Apply corrections with relaxation factor
u = u + w*du;
v = v + w*dv;
