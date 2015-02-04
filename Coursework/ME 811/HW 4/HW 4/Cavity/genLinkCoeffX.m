function [A] = genLinkCoeffX(BCu,BCv,u,v,p,rho,mu)
% Calculate the link coefficients for the x-momentum equation.

M = size(u,1);
N = size(u,2)-1;
dy=0.01/M; dx=0.01/N;

r = dy/dx;
ri = dx/dy;

%% Calculate intercell advection
cwp = @(i,j) rho*( abs(u(i,j-1)+u(i,j)) + (u(i,j-1)+u(i,j)) )/4;
cwm = @(i,j) rho*( abs(u(i,j-1)+u(i,j)) - (u(i,j-1)+u(i,j)) )/4;
cep = @(i,j) rho*( abs(u(i,j)+u(i,j+1)) + (u(i,j)+u(i,j+1)) )/4;
cem = @(i,j) rho*( abs(u(i,j)+u(i,j+1)) - (u(i,j)+u(i,j+1)) )/4;

csp = @(i,j) rho*( abs(v(i,j-1)+v(i,j)) + (v(i,j-1)+v(i,j)) )/4;
csm = @(i,j) rho*( abs(v(i,j-1)+v(i,j)) - (v(i,j-1)+v(i,j)) )/4;
cnp = @(i,j) rho*( abs(v(i+1,j-1)+v(i+1,j)) + (v(i+1,j-1)+v(i+1,j)) )/4;
cnm = @(i,j) rho*( abs(v(i+1,j-1)+v(i+1,j)) - (v(i+1,j-1)+v(i+1,j)) )/4;

%% Initialize matrices
A.s = zeros(M,N+1); 
A.w = zeros(M,N+1); 
A.o = zeros(M,N+1);
A.e = zeros(M,N+1);
A.n = zeros(M,N+1);
A.p = zeros(M,N+1);

%% Calculate inner link coefficients
i=2:M-1; j=2:N;
A.s(i,j) = -csp(i,j)*dx - mu*r;
A.w(i,j) = -cwp(i,j)*dy - mu*ri;
A.o(i,j) = cep(i,j)*dy + cwm(i,j)*dy + cnp(i,j)*dx + csm(i,j)*dx + mu*(2*r+2*ri);
A.e(i,j) = -cem(i,j)*dy - mu*ri;
A.n(i,j) = -cnm(i,j)*dx - mu*r;
A.p(i,j) = ( p(i,j-1) - p(i,j) )*dy;

%% Apply boundary conditions
% Left
i=1:M; j=1;
A.o(i,j) = 1;
A.p(i,j) = BCu.l(i);

% Right
i=1:M; j=N+1;
A.o(i,j) = 1;
A.p(i,j) = BCu.r(i);

% Bottom
i=1; j=2:N;
A.s(i,j) = 0;
A.w(i,j) = -cwp(i,j)*dy - mu*ri;
A.o(i,j) = cep(i,j)*dy + cwm(i,j)*dy + cnp(i,j)*dx + mu*(2*r+4*ri);
A.e(i,j) = -cem(i,j)*dy - mu*ri;
A.n(i,j) = -cnm(i,j)*dx - 4/3*mu*r;
A.p(i,j) = ( p(i,j-1) - p(i,j) )*dy + (rho*BCv.b(j)*dx + 8/3*mu*ri).*BCu.b(j);

% Top
i=M; j=2:N;
A.s(i,j) = -csp(i,j)*dx - 4/3*mu*r;
A.w(i,j) = -cwp(i,j)*dy - mu*ri;
A.o(i,j) = cep(i,j)*dy + cwm(i,j)*dy + csm(i,j)*dx + mu*(2*r+4*ri);
A.e(i,j) = -cem(i,j)*dy - mu*ri;
A.n(i,j) = 0;
A.p(i,j) = ( p(i,j-1) - p(i,j) )*dy + (-rho*BCv.t(j)*dx + 8/3*mu*ri).*BCu.t(j);