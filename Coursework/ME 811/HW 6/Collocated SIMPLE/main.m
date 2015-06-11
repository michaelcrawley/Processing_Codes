% Solves a pipe flow with reverse step using the SIMPLEC algorithm on a
% co-located mesh. Prescribeable boundary conditions are inlet u-velocity
% profile and outlet pressure profile.

clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = 0.02;						% Height of the pipe, m
L = 0.10;						% Length of the pipe, m

B = 5;							% Block size
M = 2*B;						% Total cells in the y-direction
N = 10*B;						% Total cells in the x-direction
m = B;							% Step cells in the y-direction
n = 2*B;						% Step cells in the x-direction

w.uv = 0.5;						% U/V-velocity relaxation factor
w.p = 0.2;						% Pressure relaxation factor
al = 0.2;						% Inertial damping factor

tol = 1e-10;					% Convergence tolerance

rho.c = ones(M,N);				% Density, kg/m3
mu.c = 2e-5*ones(M,N);			% Viscosity, kg/m-s

uin = 0.2;						% 
bc.in.u = uin*ones(M-m,1);		% Inlet u-velocity profile
bc.out.p = zeros(M,1);			% Outlet pressure profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=H/M; dx=L/N;

% Initialize matrices / initial guess
u.c = uin*ones(M,N);	u.f = uin*ones(M,N+1);
v.c = zeros(M,N);		v.f = zeros(M+1,N);
p.c = zeros(M,N);		p.f.x = zeros(M,N+1);		p.f.y = zeros(M+1,N);

% Interpolate physical values
rho.f.x = dwimX(rho.c,M,N,m,n);
rho.f.y = dwimY(rho.c,M,N,m,n);
mu.f.x = dwimX(mu.c,M,N,m,n);
mu.f.y = dwimY(mu.c,M,N,m,n);

% Prepare for first loop iteration
u.f = dwimX(u.c,M,N,m,n);
v.f = dwimY(v.c,M,N,m,n);

% Enforce boundary conditions
u.f(m+1:M,1) = bc.in.u;
p.f.x(:,N+1) = bc.out.p;
u.c(1:m,1:n) = 0;		u.f(1:m,1:n+1) = 0;
v.c(1:m,1:n) = 0;		v.f(1:m+1,1:n) = 0;

% Begin outer loop
tic;
for k=1:2000
	clc;
	fprintf('Iteration: %g\n',k);
	fprintf('Time: %g\n',toc);
	
	% Interpolate cell face pressures
	p.f.x = dwimX(p.c,M,N,m,n);
	p.f.y = dwimY(p.c,M,N,m,n);
	
	% Generate link-coefficients for x- and y-momentum equations
	links.mom = genLinkCoeffMom(bc,u,v,p,rho.f,mu.f,M,N,m,n,dy,dx);

	% Solve the x- and y-momentum equations
	[u.c R2.x(k)] = solveX(links.mom,u.c,al,4);
 	[v.c R2.y(k)] = solveY(links.mom,v.c,al,4);

	% Interpolate cell face velocities
	u.f = pwimX(bc,u.c,p.c,links.mom.o,M,N,m,n,dy);
	v.f = pwimY(v.c,p.c,links.mom.o,M,N,m,n,dx);
	
	% Calculate the mass-imbalance
	[Mimb R2.m(k)] = calcMimb(u.f,v.f,rho.f,dy,dx);
	
	% Generate link-coefficients for pressure correction equation
	links.p = genLinkCoeffP(links.mom.o,Mimb,rho.f,M,N,m,n,dy,dx);

	% Solve the pressure correction equation
	[dp R2.p(k)] = solveP(links.p,40);
	
	% Correct cell center values
	[u.c v.c] = correctUV(u.c,v.c,dp,links.mom.o,w.uv,M,N,m,n,dy,dx);
	p.c = p.c + w.p*dp;
	
	% Correct cell face velocities
	[u.f v.f] = correctUVf(u.f,v.f,dp,links.mom.o,w.uv,M,N,m,n,dy,dx);
	
	% Check for convergence
	if ( R2.x(k)<=tol && R2.y(k)<=tol && R2.p(k)<=tol ), break; end
end

fprintf('\nDone.\n');
fprintf('Completion time: %g\n',toc);