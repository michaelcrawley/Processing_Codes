clear; %close all; clc;

M = 20;
N = 20;

% Declare u-velocity boundary conditions
bc.u.t = zeros(1,N+1);				% Top
bc.u.r = sin(linspace(0,pi,M))';	% Right
bc.u.b = zeros(1,N+1);				% Bottom
bc.u.l = sin(linspace(0,pi,M))';	% Left

% Declare v-velocity boundary conditions
bc.v.t = zeros(1,N);			% Top
bc.v.r = zeros(M+1,1);			% Right
bc.v.b = zeros(1,N);			% Bottom
bc.v.l = zeros(M+1,1);			% Left

% Initialize matrices / initial guess
u = zeros(M,N+1);
v = zeros(M+1,N);
p = repmat( linspace(7,0,N) , [M 1]);

% Enforce BCs directly where possible
u(:,1) = bc.u.l;	u(:,N+1) = bc.u.r;
v(1,:) = bc.v.b;	v(M+1,:) = bc.v.t;

for k=1:1000
	A = genLinkCoeffX(bc.u,bc.v,u,v,p,1,1);

	[u R2(k)] = solveXY(A,u,0.2,4);

	figure(1); pcolor(u);
	figure(2); semilogy(R2);
	drawnow;
end