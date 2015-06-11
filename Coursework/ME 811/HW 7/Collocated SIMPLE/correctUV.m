function [u v] = correctUV(u,v,dp,Ao,w,M,N,m,n,dy,dx)

%% Initialize matrices
du = zeros(M,N);
dv = zeros(M,N);

%% Calculate u-velocity correction
% Interior cells
i=m+1:M; j=2:n+1;
du(i,j) = ( dp(i,j-1)-dp(i,j+1) )*dy./(2*Ao(i,j));
i=1:M; j=n+2:N-1;
du(i,j) = ( dp(i,j-1)-dp(i,j+1) )*dy./(2*Ao(i,j));

% Inlet
i=m+1:M; j=1;
du(i,j) = ( dp(i,j)-dp(i,j+1) )*dy./(2*Ao(i,j));

% Step face
i=1:m; j=n+1;
du(i,j) = ( dp(i,j)-dp(i,j+1) )*dy./(2*Ao(i,j));

% Outlet
i=1:M; j=N;
du(i,j) = ( dp(i,j-1)-dp(i,j) )*dy./(2*Ao(i,j));

%% Calculate v-velocity correction
% Interior cells
i=m+2:M-1; j=1:n;
dv(i,j) = ( dp(i-1,j)-dp(i+1,j) )*dx./(2*Ao(i,j));
i=2:M-1; j=n+1:N;
dv(i,j) = ( dp(i-1,j)-dp(i+1,j) )*dx./(2*Ao(i,j));

% Inlet floor
i=m+1; j=1:n;
dv(i,j) = ( dp(i,j)-dp(i+1,j) )*dx./(2*Ao(i,j));

% Main pipe floor
i=1; j=n+1:N;
dv(i,j) = ( dp(i,j)-dp(i+1,j) )*dx./(2*Ao(i,j));

% Ceiling
i=M; j=1:N;
dv(i,j) = ( dp(i-1,j)-dp(i,j) )*dx./(2*Ao(i,j));

%% Apply corrections with relaxation factor
u = u + w*du;
v = v + w*dv;
