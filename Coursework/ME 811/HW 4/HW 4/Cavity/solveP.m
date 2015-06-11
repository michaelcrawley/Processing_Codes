function [dp R2] = solveP(A,iter)
% Solve the pressure correction equation using the ADI method.

[M N] = size(A.o);

% Initial guess
null = zeros(M,N);

% Solve
dp = ADI( A.s , A.w , A.o , A.e , A.n , A.p , null , -inf , iter );

% Calculate the residual
Ri = A.s .* [zeros(1,N); dp(1:M-1,:)] ...		% phi(i,j-1)
	+ A.w .* [zeros(M,1), dp(:,1:N-1)] ...		% phi(i-1,j)
	+ A.o .* dp ...								% phi(i,j)
	+ A.e .* [dp(:,2:N), zeros(M,1)] ...		% phi(i+1,j)
	+ A.n .* [dp(2:M,:); zeros(1,N)] ...		% phi(i,j+1)
	- A.p;										% S(i,j)
R2 = sqrt(sum(sum( Ri.^2 )));