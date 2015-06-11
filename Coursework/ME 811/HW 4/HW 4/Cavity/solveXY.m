function [phi R2] = solveXY(A,phi,al,iter)
% Solve the x/y-momentum equation in correction form using the ADI method.

[M N] = size(A.o);

% Calculate the residual
Ri = A.s .* [zeros(1,N); phi(1:M-1,:)] ...		% phi(i,j-1)
	+ A.w .* [zeros(M,1), phi(:,1:N-1)] ...		% phi(i-1,j)
	+ A.o.* phi ...								% phi(i,j)
	+ A.e .* [phi(:,2:N), zeros(M,1)] ...		% phi(i+1,j)
	+ A.n .* [phi(2:M,:); zeros(1,N)] ...		% phi(i,j+1)
	- A.p;										% S(i,j)
R2 = sqrt(sum(sum( Ri.^2 )));

% Initial guess
null = zeros(M,N);

% Solve
dphi = ADI( A.s , A.w , (1+al)*A.o , A.e , A.n , -Ri , null , -inf , iter );

% Update
phi = phi + dphi;