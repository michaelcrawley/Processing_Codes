% Solves a pipe flow with reverse step using the SIMPLEC algorithm on a
% co-located mesh. Prescribeable boundary conditions are inlet u-velocity
% profile and outlet pressure profile.

clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = 0.02;						% Height of the pipe, m
L = 0.10;						% Length of the pipe, m

B = 20;							% Block size
M = 2*B;						% Total cells in the y-direction
N = 10*B;						% Total cells in the x-direction
m = B;							% Step cells in the y-direction
n = 2*B;						% Step cells in the x-direction

w.uv = 0.5;						% U/V-velocity relaxation factor
w.p = 0.2;						% Pressure relaxation factor
al = 0.2;						% Inertial damping factor

tol = 1e-10;					% Convergence tolerance

mu.c = 2e-5*ones(M,N);			% Viscosity, kg/m-s
cp.c = 1012*ones(M,N);			% Heat capacity, J/kg-K
k.c = 0.026*ones(M,N);			% Thermal conductivity, J/m-K
R = 287;						% Specific gas constant, ???
Tref = 300;						% Reference temperature, K
pref = 101325;					% Reference pressure, Pa

uin = 0.1;						% 
bc.in.u = uin*ones(M-m,1);		% Inlet u-velocity profile, m/s
bc.out.p = zeros(M,1);			% Outlet pressure profile, ??
bc.stop.T = 100;				% Step top temperature, K
bc.sface.T = 100;				% Step face temperature, K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=H/M; dx=L/N;

% Initialize matrices / initial guess
u.c = uin*ones(M,N);	u.f = uin*ones(M,N+1);
v.c = zeros(M,N);		v.f = zeros(M+1,N);
p.c = zeros(M,N);		p.f.x = zeros(M,N+1);		p.f.y = zeros(M+1,N);
T.c = zeros(M,N);

% Interpolate physical properties
mu.f.x = dwimX(mu.c,M,N,m,n);		mu.f.y = dwimY(mu.c,M,N,m,n);
cp.f.x = dwimX(cp.c,M,N,m,n);		cp.f.y = dwimY(cp.c,M,N,m,n);
k.f.x = dwimX(k.c,M,N,m,n);			k.f.y = dwimY(k.c,M,N,m,n);

% Prepare for first loop iteration
u.f = dwimX(u.c,M,N,m,n);
v.f = dwimY(v.c,M,N,m,n);
rho.c = calcRho(p.c,T.c,pref,Tref,R);

% Enforce boundary conditions
u.f(m+1:M,1) = bc.in.u;				p.f.x(:,N+1) = bc.out.p;
u.c(1:m,1:n) = 0;					u.f(1:m,1:n+1) = 0;
v.c(1:m,1:n) = 0;					v.f(1:m+1,1:n) = 0;

% Begin outer loop
tic;
for kk=1:2000
	clc;
	fprintf('Iteration: %g\n',kk);
	fprintf('Time: %g\n',toc);
	
	% Interpolate cell face values
	p.f.x = dwimX(p.c,M,N,m,n);			p.f.y = dwimY(p.c,M,N,m,n);
	rho.f.x = dwimX(rho.c,M,N,m,n);		rho.f.y = dwimY(rho.c,M,N,m,n);
	
	% Generate link-coefficients for x- and y-momentum equations
	links.mom = genLinkCoeffMom(bc,u,v,p,rho.f,mu.f,M,N,m,n,dy,dx);

	% Solve the x- and y-momentum equations
	[u.c R2.x(kk)] = solveX(links.mom,u.c,al,4);
 	[v.c R2.y(kk)] = solveY(links.mom,v.c,al,4);

	% Interpolate cell face velocities
	u.f = pwimX(bc,u.c,p.c,links.mom.o,M,N,m,n,dy);
	v.f = pwimY(v.c,p.c,links.mom.o,M,N,m,n,dx);
	
	% Calculate the mass-imbalance
	[Mimb R2.m(kk)] = calcMimb(u.f,v.f,rho.f,dy,dx);
	
	% Generate link-coefficients for pressure correction equation
	links.p = genLinkCoeffP(links.mom.o,Mimb,rho.f,M,N,m,n,dy,dx);

	% Solve the pressure correction equation
	[dp R2.p(kk)] = solveP(links.p,40);
	
	% Correct cell center values
	[u.c v.c] = correctUV(u.c,v.c,dp,links.mom.o,w.uv,M,N,m,n,dy,dx);
	p.c = p.c + w.p*dp;
	
	% Correct cell face velocities
	[u.f v.f] = correctUVf(u.f,v.f,dp,links.mom.o,w.uv,M,N,m,n,dy,dx);
	
	% Generate link-coefficients for energy equation
	links.T = genLinkCoeffT(bc,u,v,p,rho.f,cp.f,k.f,M,N,m,n,dy,dx);
	
	% Solve the energy equation
	[T.c R2.T(kk)] = solveT(links.T,T.c,al,40);
	
	% Correct cell center density
	rho.c = calcRho(p.c,T.c,pref,Tref,R);
	
	% Check for convergence
	if ( R2.x(kk)<=tol && R2.y(kk)<=tol && R2.p(kk)<=tol && R2.T(kk)<=tol ), break; end
end

fprintf('\nDone.\n');
fprintf('Completion time: %g\n',toc);
