function [A] = genLinkCoeffP(Ao,Mimb,rhof,M,N,m,n,dy,dx)
% Calculate the link coefficients for the pressure correction equation.

%% Initialize matrices
A.s = zeros(M,N); 
A.w = zeros(M,N); 
A.o = zeros(M,N);
A.e = zeros(M,N);
A.n = zeros(M,N);
A.p = zeros(M,N);

%% Calculate interior link coefficients
% Inlet pipe
i=m+2:M-1; j=2:n+2;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Main pipe
i=2:M-1; j=n+2:N-1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Awkward corner
i=m+1; j=n+1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

%% Calculate boundary link coefficients
% Inlet
i=m+2:M-1; j=1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = 0;
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Outlet
i=2:M-1; j=N;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = 0;
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j)) - rhof.x(i,j+1)*dy^2.*(1./Ao(i,j));
A.p(i,j) = Mimb(i,j);

% Top
i=M; j=2:N-1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom (inlet)
i=m+1; j=2:n;
A.s(i,j) = 0;
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom (main pipe)
i=1; j=n+2:N-1;
A.s(i,j) = 0;
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Step face
i=2:m; j=n+1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = 0;
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

%% Calculate boundary link coefficients for corner cells
% Top left
i=M; j=1;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = 0;
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Top right
i=M; j=N;
A.s(i,j) = rhof.y(i,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i-1,j));
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = 0;
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j)) - rhof.x(i,j+1)*dy^2.*(1./Ao(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom left
i=m+1; j=1;
A.s(i,j) = 0;
A.w(i,j) = 0;
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom right
i=1; j=N;
A.s(i,j) = 0;
A.w(i,j) = rhof.x(i,j)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j-1));
A.e(i,j) = 0;
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j)) - rhof.x(i,j+1)*dy^2.*(1./Ao(i,j));
A.p(i,j) = Mimb(i,j);

% Step bottom
i=1; j=n+1;
A.s(i,j) = 0;
A.w(i,j) = 0;
A.e(i,j) = rhof.x(i,j+1)*dy^2/2.*(1./Ao(i,j)+1./Ao(i,j+1));
A.n(i,j) = rhof.y(i+1,j)*dx^2/2.*(1./Ao(i,j)+1./Ao(i+1,j));
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

%% Populate step interior to ensure well defined system
i=1:m; j=1:n;
A.o(i,j) = 1;