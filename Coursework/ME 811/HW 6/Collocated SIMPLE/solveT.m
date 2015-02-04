function [T R2] = solveT(A,T,al,iter)
% Solve the energy equation in correction form using the ADI method.

[M N] = size(A.o);

% Calculate the residual
Ri = A.s .* [zeros(1,N); T(1:M-1,:)] ...		% T(i,j-1)
	+ A.w .* [zeros(M,1), T(:,1:N-1)] ...		% T(i-1,j)
	+ A.o.* T ...								% T(i,j)
	+ A.e .* [T(:,2:N), zeros(M,1)] ...			% T(i+1,j)
	+ A.n .* [T(2:M,:); zeros(1,N)] ...			% T(i,j+1)
	- A.p;										% S(i,j)
R2 = sqrt(sum(sum( Ri.^2 )));

% Initial guess
null = zeros(M,N);

% Solve
dT = ADI( A.s , A.w , (1+al)*A.o , A.e , A.n , -Ri , null , -inf , iter );

% Update
T = T + dT;