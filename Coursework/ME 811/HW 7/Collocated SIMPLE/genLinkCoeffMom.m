function [A] = genLinkCoeffMom(bc,u,v,p,rhof,muf,M,N,m,n,dy,dx)

ra = dy/dx;
ar = dx/dy;

%% Calculate intercell advection
cwp = @(i,j) rhof.x(i,j).*( abs(u.f(i,j))+u.f(i,j) )/2;
cwm = @(i,j) rhof.x(i,j).*( abs(u.f(i,j))-u.f(i,j) )/2;
cep = @(i,j) rhof.x(i,j+1).*( abs(u.f(i,j+1))+u.f(i,j+1) )/2;
cem = @(i,j) rhof.x(i,j+1).*( abs(u.f(i,j+1))-u.f(i,j+1) )/2;

csp = @(i,j) rhof.y(i,j).*( abs(v.f(i,j))+v.f(i,j) )/2;
csm = @(i,j) rhof.y(i,j).*( abs(v.f(i,j))-v.f(i,j) )/2;
cnp = @(i,j) rhof.y(i+1,j).*( abs(v.f(i+1,j))+v.f(i+1,j) )/2;
cnm = @(i,j) rhof.y(i+1,j).*( abs(v.f(i+1,j))-v.f(i+1,j) )/2;

%% Initialize matrices
A.s = zeros(M,N);
A.w = zeros(M,N);
A.o = zeros(M,N);
A.e = zeros(M,N);
A.n = zeros(M,N);
A.p.x = zeros(M,N);
A.p.y = zeros(M,N);

%% Calculate interior link coefficients
% Inlet pipe
i=m+2:M-1; j=2:n+1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

% Main pipe
i=2:M-1; j=n+2:N-1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

% Awkward cell at step corner
i=m+1; j=n+1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

%% Calculate boundary link coefficients
% Inlet
i=m+2:M-1; j=1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.o(i,j) = cep(i,j)*dy + muf.x(i,j+1)*ar + 3*muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar - muf.x(i,j)*ar/3;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
A.p.x(i,j) = (p.c(i,j)-p.f.x(i,j+1))*dy ...
	+ (rhof.x(i,j).*bc.in.u(2:end-1)*dy + 8/3*muf.x(i,j)*ra).*bc.in.u(2:end-1);
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

% Outlet
i=2:M-1; j=N;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (rhof.x(i,j+1).*u.c(i,j)+cwm(i,j))*dy + muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
% A.p.x(i,j) = -bc.out.p(i)*dy;
A.p.x(i,j) = (p.f.x(i,j)-bc.out.p(i))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

% Top
i=M; j=2:N-1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra - muf.y(i+1,j)*ra/3;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ csm(i,j)*dx + 3*muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.c(i,j))*dx;

% Bottom (inlet)
i=m+1; j=2:n;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ cnp(i,j)*dx + muf.y(i+1,j)*ra + 3*muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra - muf.y(i,j)*ra/3;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.c(i,j)-p.f.y(i+1,j))*dx;

% Bottom (main pipe)
i=1; j=n+2:N-1;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (cep(i,j)+cwm(i,j))*dy + muf.x(i,j+1)*ar + muf.x(i,j)*ar ...
	+ cnp(i,j)*dx + muf.y(i+1,j)*ra + 3*muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra - muf.y(i,j)*ra/3;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.c(i,j)-p.f.y(i+1,j))*dx;

% Step face
i=2:m; j=n+1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra;
A.o(i,j) = cep(i,j)*dy + muf.x(i,j+1)*ar + 3*muf.x(i,j)*ar ...
	+ (cnp(i,j)+csm(i,j))*dx + muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar - muf.x(i,j)*ar/3;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra;
A.p.x(i,j) = (p.c(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.f.y(i+1,j))*dx;

%% Calculate boundary link coefficients for corner cells
% Top left
i=M; j=1;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra - muf.y(i+1,j)*ra/3;
A.o(i,j) = cep(i,j)*dy + muf.x(i,j+1)*ar + 3*muf.x(i,j)*ar ...
	+ csm(i,j)*dx + 3*muf.y(i+1,j)*ra + muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar - muf.x(i,j)*ar/3;
A.p.x(i,j) = (p.f.x(i,j)-p.f.x(i,j+1))*dy ...
	+ (rhof.x(i,j).*bc.in.u(end)*dy + 8/3*muf.x(i,j)*ra)*bc.in.u(end);
A.p.y(i,j) = 0;

% Top right
i=M; j=N;
A.s(i,j) = -csp(i,j)*dx - muf.y(i,j)*ra - muf.y(i+1,j)*ra/3;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (rhof.x(i,j+1).*u.c(i,j)+cwm(i,j))*dy + muf.x(i,j)*ar ...
	+ csm(i,j)*dx + 3*muf.y(i+1,j)*ra + muf.y(i,j)*ra;
% A.p.x(i,j) = -bc.out.p(i)*dy;
A.p.x(i,j) = (p.f.x(i,j)-bc.out.p(i))*dy;
A.p.y(i,j) = (p.f.y(i,j)-p.c(i,j))*dx;

% Bottom left
i=m+1; j=1;
A.o(i,j) = cep(i,j)*dy + muf.x(i,j+1)*ar + 3*muf.x(i,j)*ar ...
	+ cnp(i,j)*dx + muf.y(i+1,j)*ra + 3*muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar - muf.x(i,j)*ar/3;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra - muf.y(i,j)*ra/3;
A.p.x(i,j) = (p.c(i,j)-p.f.x(i,j+1))*dy ...
	+ (rhof.x(i,j).*bc.in.u(1)*dy +8/3*muf.x(i,j)*ra)*bc.in.u(1);
A.p.y(i,j) = (p.c(i,j)-p.f.y(i+1,j))*dx;

% Bottom right
i=1; j=N;
A.w(i,j) = -cwp(i,j)*dy - muf.x(i,j)*ar;
A.o(i,j) = (rhof.x(i,j+1).*u.c(i,j)+cwm(i,j))*dy + muf.x(i,j)*ar ...
	+ cnp(i,j)*dx + muf.y(i+1,j)*ra + 3*muf.y(i,j)*ra;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra - muf.y(i,j)*ra/3;
% A.p.x(i,j) = -bc.out.p(i)*dy;
A.p.x(i,j) = (p.f.x(i,j)-bc.out.p(i))*dy;
A.p.y(i,j) = (p.c(i,j)-p.f.y(i+1,j))*dx;

% Step bottom
i=1; j=n+1;
A.o(i,j) = cep(i,j)*dy + muf.x(i,j+1)*ar + 3*muf.x(i,j)*ar ...
	+ cnp(i,j)*dx + muf.y(i+1,j)*ra + 3*muf.y(i,j)*ra;
A.e(i,j) = -cem(i,j)*dy - muf.x(i,j+1)*ar - muf.x(i,j)*ar/3;
A.n(i,j) = -cnm(i,j)*dx - muf.y(i+1,j)*ra - muf.y(i,j)*ra/3;
A.p.x(i,j) = (p.c(i,j)-p.f.x(i,j+1))*dy;
A.p.y(i,j) = (p.c(i,j)-p.f.y(i+1,j))*dx;

%% Populate step interior to ensure well defined system
i=1:m; j=1:n;
A.o(i,j) = 1;