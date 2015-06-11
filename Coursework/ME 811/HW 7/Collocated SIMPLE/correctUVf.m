function [uf vf] = correctUVf(uf,vf,dp,Ao,w,M,N,m,n,dy,dx)

%% Initialize matrices
duf = zeros(M,N+1);
dvf = zeros(M+1,N);

%% Calculate corrections
i=m+1:M; j=2:n+1;
duf(i,j) = dy/2*(1./Ao(i,j)+1./Ao(i,j-1)).*(dp(i,j-1)-dp(i,j));
i=1:M; j=n+2:N;
duf(i,j) = dy/2*(1./Ao(i,j)+1./Ao(i,j-1)).*(dp(i,j-1)-dp(i,j));

i=m+2:M; j=1:n;
dvf(i,j) = dx/2*(1./Ao(i,j)+1./Ao(i-1,j)).*(dp(i-1,j)-dp(i,j));
i=2:M; j=n+1:N;
dvf(i,j) = dx/2*(1./Ao(i,j)+1./Ao(i-1,j)).*(dp(i-1,j)-dp(i,j));

%% Apply corrections with relaxation factor
uf = uf + w*duf;
vf = vf + w*dvf;