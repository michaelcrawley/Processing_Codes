function [A] = genLinkCoeffP(Aox,Aoy,Mimb,rho)
% Calculate the link coefficients for the pressure correction equation.

[M N] = size(Mimb);
dy=0.01/M; dx=0.01/N;

%% Initialize matrices
A.s = zeros(M,N); 
A.w = zeros(M,N); 
A.o = zeros(M,N);
A.e = zeros(M,N);
A.n = zeros(M,N);
A.p = zeros(M,N);

%% Calculate inner link coefficients
i=2:M-1; j=2:N-1;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

%% Apply boundary conditions
% Top
i=M; j=2:N-1;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom
i=1; j=2:N-1;
A.s(i,j) = 0;
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Left
i=2:M-1; j=1;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = 0;
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Right
i=2:M-1; j=N;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = 0;
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Top-left
i=M; j=1;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = 0;
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Top-right
i=M; j=N;
A.s(i,j) = rho*dx^2./Aoy(i,j);
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = 0;
A.n(i,j) = 0;
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom-left
i=1; j=1;
A.s(i,j) = 0;
A.w(i,j) = 0;
A.e(i,j) = rho*dy^2./Aox(i,j+1);
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);

% Bottom-right
i=1; j=N;
A.s(i,j) = 0;
A.w(i,j) = rho*dy^2./Aox(i,j);
A.e(i,j) = 0;
A.n(i,j) = rho*dx^2./Aoy(i+1,j);
A.o(i,j) = -(A.s(i,j) + A.w(i,j) + A.e(i,j) + A.n(i,j));
A.p(i,j) = Mimb(i,j);