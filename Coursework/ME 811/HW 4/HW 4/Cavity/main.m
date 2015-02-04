clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 40;							% Cells (rows) in the y-direction
N = 40;							% Cells (columns) in the x-direction

lax.uv = 0.5;					% U/V-velocity relaxation factor
lax.p = 0.2;					% Pressure relaxation factor
al = 0.2;						% Inertial damping factor

tol = 1e-10;					% Convergence tolerance

rho = 1000;						% Density, kg/m3
mu = 1e-3;						% Viscosity, kg/m-s

ulid = 0.01;					% Lid velocity, 0.01 ~ Re 100

% Declare u-velocity boundary conditions
bc.u.t = ulid*ones(1,N+1);		% Top
bc.u.r = zeros(M,1);			% Right
bc.u.b = zeros(1,N+1);			% Bottom
bc.u.l = zeros(M,1);			% Left

% Declare v-velocity boundary conditions
bc.v.t = zeros(1,N);			% Top
bc.v.r = zeros(M+1,1);			% Right
bc.v.b = zeros(1,N);			% Bottom
bc.v.l = zeros(M+1,1);			% Left

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=0.01/M; dx=0.01/N;

% Initialize matrices / initial guess
u = zeros(M,N+1);
v = zeros(M+1,N);
p = zeros(M,N);

% Enforce BCs where direct overlay exists
u(:,1) = bc.u.l;	u(:,N+1) = bc.u.r;
v(1,:) = bc.v.b;	v(M+1,:) = bc.v.t;

for k=1:2000
	clc;
	fprintf('Iteration: %g\n',k);
	
	fprintf('X+Y momentum link gen... ');
	% Generate link-coefficients for x- and y-momentum equations
	links.x = genLinkCoeffX(bc.u,bc.v,u,v,p,rho,mu);
 	links.y = genLinkCoeffY(bc.u,bc.v,u,v,p,rho,mu);
	fprintf('done.\n');
	
	fprintf('X+Y momentum solution... ');
	% Solve the x- and y-momentum equations
	[u R2.x(k)] = solveXY(links.x,u,al,4);
 	[v R2.y(k)] = solveXY(links.y,v,al,4);
	fprintf('done.\n');
	
	fprintf('Mass imbalance calculation... ');
	% Calculate the mass-imbalance
	[Mimb R2.m(k)] = calcMimb(u,v,rho);
	fprintf('done.\n');
	
	fprintf('P correction link gen... ');
	% Generate link-coefficients for pressure correction equation
	links.p = genLinkCoeffP(links.x.o,links.y.o,Mimb,rho);
	fprintf('done.\n');
	
	fprintf('P correction solution... ');
	% Solve the pressure correction equation
	[dp R2.p(k)] = solveP(links.p,40);
	fprintf('done.\n');
	
	fprintf('U, V, P update... ');
	% Update velocity and pressure fields
	[u v] = correctUV(u,v,dp,links.x.o,links.y.o,lax.uv);
	p = p + lax.p*dp;
	fprintf('done.\n');
		
	% Check for convergence
	if ( R2.x(k)<=tol && R2.y(k)<=tol && R2.p(k)<=tol ), break; end
end

% Interpolate velocity fields to cell centers
uc = ( u(:,1:N) + u(:,2:N+1) )/2;
vc = ( v(1:M,:) + v(2:M+1,:) )/2;
[x y] = meshgrid( dx/2:dx:0.01-dx/2 , dy/2:dy:0.01-dy/2 );

% % % Save outputs
% % save('n40re100.mat');
% % 
% % fprintf('\nDone.\n');